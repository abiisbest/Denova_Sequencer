import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="De Novo Professional Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "type": "Gram-Negative"},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "type": "Gram-Positive"},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "type": "Gram-Positive"},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "type": "Eukaryote"},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "type": "Acid-Fast"}
}

def remove_adapters(reads, adapter_seq, min_keep_len):
    cleaned_reads = []
    for read in reads:
        if adapter_seq and adapter_seq in read:
            read = read.split(adapter_seq)[0]
        if len(read) >= min_keep_len:
            cleaned_reads.append(read)
    return cleaned_reads

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def annotate_gene(seq, species_gc):
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    length = len(seq)
    
    if length > 1200:
        return "Polymerase/Replication Factor"
    elif abs(gc - species_gc) > 5:
        return "Horizontal Gene Transfer (HGT) Candidate"
    elif length > 600:
        return "Metabolic Enzyme"
    else:
        return "Hypothetical Protein"

def find_all_orfs(sequence, species_gc, min_len=300):
    found_genes = []
    pattern = re.compile(r'(ATG(?:...){%d,1000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    for strand_name, strand_sign in [("Forward", "+"), ("Reverse", "-")]:
        dna = sequence if strand_name == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                start_pos = match.start() + frame
                
                found_genes.append({
                    "Locus_Tag": f"GENE_{len(found_genes)+1:04d}",
                    "Strand": strand_sign,
                    "Start": int(start_pos),
                    "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)),
                    "GC_%": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Function": annotate_gene(gene_seq, species_gc),
                    "Sequence": gene_seq
                })
    return found_genes

st.title("🧬 De Novo: Auto-Identifying Genomic Pipeline")
st.markdown("---")

st.sidebar.header("⚙️ Pipeline Settings")
trim_adapters = st.sidebar.checkbox("Enable Adapter Trimming", value=True)
adapter_sequence = st.sidebar.text_input("Adapter Sequence", "AGATCGGAAGAG")
min_len_filter = st.sidebar.slider("Minimum Read Length (bp)", 0, 100, 15)

uploaded_file = st.file_uploader("Upload Raw FASTQ/GZ Data", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        raw_reads = [line.strip() for line in data.splitlines()[1::4] if line.strip()]
        
        if not raw_reads:
            st.error("No valid DNA reads found.")
            st.stop()

        sample_gc = round(("".join(raw_reads[:100]).count('G') + "".join(raw_reads[:100]).count('C')) / len("".join(raw_reads[:100])) * 100, 2)
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.info(f"Predicted Species: **{closest_species}** (GC: {sample_gc}%)")
        ref = SPECIES_LIBRARY[closest_species]

        if st.button("🚀 Run Full Annotation"):
            trimmed_reads = remove_adapters(raw_reads, adapter_sequence, min_len_filter) if trim_adapters else raw_reads
            
            if not trimmed_reads:
                st.error("Trimming removed all reads. Lower the 'Minimum Read Length' in sidebar.")
                st.stop()

            full_genome = "NNNNN".join(trimmed_reads[:200])
            all_orfs = find_all_orfs(full_genome, ref['ref_gc'])
            genes_df = pd.DataFrame(all_orfs).sort_values('Start').drop_duplicates(subset=['Start'])

            tab1, tab2 = st.tabs(["📊 Quality Control", "🧬 Functional Annotation"])

            with tab1:
                col1, col2 = st.columns(2)
                col1.metric("Reads Processed", len(trimmed_reads))
                col2.metric("Morphology", ref['type'])
                if "Gram-Positive" in ref['type']:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/8/86/Gram_positive_cell_wall.svg/300px-Gram_positive_cell_wall.svg.png")
                else:
                    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/d/d3/Gram_negative_cell_wall.svg/300px-Gram_negative_cell_wall.svg.png")

            with tab2:
                st.subheader("📝 Gene Feature Table")
                st.dataframe(genes_df.drop(columns=['Sequence']), use_container_width=True)
                
                gff = "##gff-version 3\n" + "\n".join([f"contig1\tDeNovo\tCDS\t{r['Start']}\t{r['End']}\t.\t{r['Strand']}\t0\tID={r['Locus_Tag']};Note={r['Function']}" for i, r in genes_df.iterrows()])
                st.download_button("💾 Download GFF3 Annotation", gff, "annotation.gff3")

    except Exception as e:
        st.error(f"Error: {e}")

import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import io

st.set_page_config(page_title="De Novo Professional Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "expected_genes": 4300, "genome_size": 4641652},
    "Staphylococcus aureus": {"ref_gc": 32.8, "expected_genes": 2800, "genome_size": 2821361},
    "Bacillus subtilis": {"ref_gc": 43.5, "expected_genes": 4200, "genome_size": 4214630},
    "Saccharomyces cerevisiae (Yeast)": {"ref_gc": 38.3, "expected_genes": 6000, "genome_size": 12157105},
    "Mycobacterium tuberculosis": {"ref_gc": 65.6, "expected_genes": 4000, "genome_size": 4411532}
}

st.title("🧬 De Novo: Auto-Identifying Genomic Pipeline")
st.markdown("---")

def format_indian_num(n):
    if n >= 100000:
        return f"{n/100000:.2f}L"
    return f"{n:,}"

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300):
    found_genes = []
    pattern = re.compile(r'(ATG(?:...){%d,1000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                start_pos = match.start() + frame
                found_genes.append({
                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Sequence": gene_seq, "Type": "Gene"
                })
    return found_genes

def parse_genomic_file(uploaded_file):
    filename = uploaded_file.name
    if filename.endswith('.gz'):
        content = gzip.decompress(uploaded_file.read()).decode("utf-8")
    else:
        content = uploaded_file.read().decode("utf-8")
    
    lines = content.splitlines()
    reads = []
    
    if filename.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        reads = [line.strip() for line in lines[1::4] if line.strip()]
    elif filename.endswith(('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz')):
        current_seq = []
        for line in lines:
            if line.startswith(">"):
                if current_seq:
                    reads.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            reads.append("".join(current_seq))
    else:
        reads = [line.strip() for line in lines if line.strip() and not line.startswith(('@', '>', '+'))]
        
    return reads

uploaded_file = st.file_uploader("Upload Genomic Data (FASTA, FASTQ, GZ)", type=["fastq", "fq", "fasta", "fa", "fna", "gz"])

if uploaded_file:
    try:
        raw_reads = parse_genomic_file(uploaded_file)
        
        if not raw_reads:
            st.error("No valid sequences found in the file.")
            st.stop()

        sample_seq = "".join(raw_reads[:500])
        sample_gc = round((sample_seq.count('G') + sample_seq.count('C')) / len(sample_seq) * 100, 2)
        
        closest_species = min(SPECIES_LIBRARY.keys(), key=lambda x: abs(SPECIES_LIBRARY[x]['ref_gc'] - sample_gc))
        
        st.subheader("🤖 Auto-Identification Result")
        st.info(f"Sample GC Content: **{sample_gc}%** | Predicted Species: **{closest_species}**")
        
        is_correct = st.radio("Confirm Species Identity:", ("Yes, proceed with prediction", "No, let me choose manually"))
        selected_species = closest_species if is_correct == "Yes, proceed with prediction" else st.selectbox("Select Species", list(SPECIES_LIBRARY.keys()))
        ref = SPECIES_LIBRARY.get(selected_species)

        if st.button("🚀 Run Full Analysis"):
            raw_len_avg = sum(len(r) for r in raw_reads) / len(raw_reads)
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60] if "fastq" in uploaded_file.name else raw_reads
            trim_len_avg = sum(len(r) for r in trimmed_reads) / len(trimmed_reads)
            
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            all_raw_orfs = find_all_orfs(full_genome)
            genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC Report", "🏗️ Reference Alignment", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Data Statistics")
                col_qc1, col_qc2 = st.columns(2)
                with col_qc1:
                    st.metric("Total Contigs/Reads", format_indian_num(len(raw_reads)))
                    st.metric("Avg Length", f"{raw_len_avg:.1f} bp")
                with col_qc2:
                    st.metric("Total Bases Processed", format_indian_num(sum(len(r) for r in trimmed_reads)))
                    st.metric("Max Sequence Length", f"{max(len(r) for r in raw_reads):,} bp")

                results_data = {
                    "Metric Parameter": ["Genomic GC Signature", "ORF Discovery Yield", "Coding Density", "Reference Conformity"],
                    "Observed Result": [f"{sample_gc}%", f"{len(genes_df)} features", f"{round((genes_df['Length'].sum()/total_len)*100, 2)}%", f"{round(100 - abs(sample_gc - ref['ref_gc']), 2)}%"]
                }
                st.table(pd.DataFrame(results_data))

            with tab2:
                st.subheader("🏗️ Comparative Alignment Analysis")
                current_gc = round((full_genome.count('G') + full_genome.count('C')) / len(full_genome) * 100, 2)
                st.metric("Sample GC %", f"{current_gc}%", f"{current_gc - ref['ref_gc']:.2f}% Dev from Ref")

                window = max(100, total_len // 100)
                p_skew, skews = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)

                fig_skew = go.Figure()
                skews_np = np.array(skews)
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np >= 0, skews_np, 0), fill='tozeroy', mode='lines', line=dict(color='#00CC96', width=0), name='Pos Skew'))
                fig_skew.add_trace(go.Scatter(x=p_skew, y=np.where(skews_np < 0, skews_np, 0), fill='tozeroy', mode='lines', line=dict(color='#EF553B', width=0), name='Neg Skew'))
                fig_skew.update_layout(title="GC Skew Analysis", template="plotly_dark", height=400)
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🧬 Annotation Performance")
                st.dataframe(genes_df.drop(columns=['Sequence', 'Type']), use_container_width=True)
                
                fasta_out = "".join([f">gene_{i}\n{r['Sequence']}\n" for i, r in genes_df.iterrows()])
                st.download_button("📝 Download FASTA Proteins", fasta_out, "predicted_genes.fasta")

    except Exception as e:
        st.error(f"Analysis Error: {e}")

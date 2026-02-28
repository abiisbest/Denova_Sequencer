import streamlit as st
import gzip
import re
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(page_title="De Nova Professional", layout="wide")

st.title("🧬 De Nova: Professional Genome Assembly & Analysis")
st.markdown("---")

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
                    "Strand": strand,
                    "Start": start_pos,
                    "End": start_pos + len(gene_seq),
                    "Length": len(gene_seq),
                    "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2)
                })
    return found_genes

uploaded_file = st.file_uploader("Upload FASTQ or GZ File", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        reads = [line.strip() for line in data.splitlines()[1::4] if len(line) > 50]

        if st.button("🚀 Run Full Genomic Pipeline"):
            full_genome = "NNNNN".join(reads[:200]) 
            total_len = len(full_genome)
            
            # --- 1. METRICS ---
            st.subheader("📊 Assembly Quality Metrics")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Length", f"{total_len} bp")
            col2.metric("Avg GC %", f"{round((full_genome.count('G')+full_genome.count('C'))/total_len*100, 2)}%")
            col3.metric("Contigs Used", len(reads[:200]))

            # --- 2. THE MAIN GRAPH: GC SKEW ---
            st.subheader("📈 GC Skew Analysis (Origin of Replication)")
            
            window = 500
            skews, positions = [], []
            for i in range(0, total_len - window, window):
                sub = full_genome[i:i+window]
                g, c = sub.count('G'), sub.count('C')
                skew = (g - c) / (g + c) if (g + c) > 0 else 0
                skews.append(skew)
                positions.append(i)
            
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(positions, skews, color='#1f77b4', label='GC Skew')
            ax.axhline(0, color='red', linestyle='--', alpha=0.5)
            ax.set_xlabel("Genome Position (bp)")
            ax.set_ylabel("Skew (G-C)/(G+C)")
            ax.grid(True, alpha=0.3)
            st.pyplot(fig)

            # --- 3. ANNOTATION TABLE ---
            st.subheader("🧬 Predicted Coding Sequences (CDS)")
            all_genes = find_all_orfs(full_genome)
            if all_genes:
                df = pd.DataFrame(all_genes).sort_values('Start')
                df = df.drop_duplicates(subset=['Start', 'Strand'], keep='first')
                st.dataframe(df, use_container_width=True)
                
                gff = "##gff-version 3\n"
                for i, row in df.iterrows():
                    s = "+" if row['Strand'] == "Forward" else "-"
                    gff += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{s}\t0\tID=gene_{i}\n"
                st.download_button("💾 Download GFF3 Annotation", gff, "annotation.gff3")

    except Exception as e:
        st.error(f"Error: {e}")

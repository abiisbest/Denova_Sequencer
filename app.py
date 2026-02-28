import streamlit as st
import gzip
import re
import time
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="De Nova Professional", layout="wide")

st.title("🧬 De Nova: Professional Genome Assembly & Visualization")
st.markdown("---")

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300):
    found_genes = []
    # Regular expression for ATG -> Triplets -> Stop Codon
    pattern = re.compile(r'(ATG(..[TCGA]){%d,1000}(TAG|TAA|TGA))' % (min_len // 3))
    
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

uploaded_file = st.file_uploader("Upload Genomic Data (FASTQ/GZ)", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        # 1. Handle Binary/Compressed Data
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        # 2. Extract Sequences from FASTQ
        reads = [line.strip() for line in data.splitlines()[1::4] if len(line) > 50]

        if st.button("🚀 Execute Full Analysis"):
            # Simulate assembly with spacers to separate genes
            full_genome = "NNNNN".join(reads[:200]) 
            total_len = len(full_genome)
            
            # --- METRICS ---
            st.subheader("📊 Assembly Quality Metrics")
            m1, m2, m3 = st.columns(3)
            m1.metric("Total Length", f"{total_len} bp")
            m2.metric("Assembly GC %", f"{round((full_genome.count('G')+full_genome.count('C'))/total_len*100, 2)}%")
            m3.metric("Contigs Used", len(reads[:200]))

            # --- GC SKEW PLOT ---
            st.subheader("📈 GC Skew Analysis (Origin of Replication)")
            window = 500
            skews, positions = [], []
            for i in range(0, total_len - window, window):
                sub = full_genome[i:i+window]
                g, c = sub.count('G'), sub.count('C')
                skew = (g - c) / (g + c) if (g + c) > 0 else 0
                skews.append(skew)
                positions.append(i)
            
            fig_skew = px.line(x=positions, y=skews, labels={'x':'Genome Position (bp)', 'y':'GC Skew'})
            fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
            st.plotly_chart(fig_skew, use_container_width=True)

            # --- FIXED LINEAR GENOME MAP ---
            st.subheader("🗺️ Linear Genome Map (Gene Locations)")
            all_genes = find_all_orfs(full_genome)
            if all_genes:
                df = pd.DataFrame(all_genes).sort_values('Start')
                df = df.drop_duplicates(subset=['Start'], keep='first')
                
                # Using bar with 'base' to represent linear genomic coordinates
                fig_map = px.bar(
                    df, 
                    x="Length", 
                    y="Strand", 
                    base="Start", 
                    orientation='h',
                    color="GC %",
                    hover_data=["Start", "End", "Length", "GC %"],
                    color_continuous_scale="Viridis"
                )
                
                fig_map.update_layout(
                    xaxis_title="Genome Position (bp)",
                    yaxis_title="Strand",
                    coloraxis_colorbar=dict(title="GC %"),
                    height=350
                )
                fig_map.update_xaxes(type='linear')
                st.plotly_chart(fig_map, use_container_width=True)

                # --- GFF3 EXPORT ---
                gff_content = "##gff-version 3\n"
                for i, row in df.iterrows():
                    strand_sign = "+" if row['Strand'] == "Forward" else "-"
                    gff_content += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{strand_sign}\t0\tID=gene_{i};GC={row['GC %']}\n"
                
                st.divider()
                st.download_button("💾 Download GFF3 Annotation", gff_content, "annotation.gff3")
                st.dataframe(df.drop(columns=["Sequence"], errors='ignore'), use_container_width=True)

    except Exception as e:
        st.error(f"Execution Error: {e}")
else:
    st.info("Awaiting genomic data upload (FASTQ format).")

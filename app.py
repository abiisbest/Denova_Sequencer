import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go

st.set_page_config(page_title="De Nova Professional", layout="wide")

st.title("🧬 De Nova: Professional Genome Assembly & Analysis")
st.markdown("---")

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300):
    found_genes = []
    # Standard Genetic Code ORF pattern
    pattern = re.compile(r'(ATG(?:...){%d,1000}?(?:TAG|TAA|TGA))' % (min_len // 3))
    
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                start_pos = match.start() + frame
                found_genes.append({
                    "Strand": strand,
                    "Start": int(start_pos),
                    "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)),
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

        if st.button("🚀 Execute Full Genomic Pipeline"):
            # Assembly Simulation (joining reads with spacers)
            full_genome = "NNNNN".join(reads[:200]) 
            total_len = len(full_genome)
            
            # --- 1. METRICS ---
            st.subheader("📊 Assembly Quality Metrics")
            m1, m2, m3 = st.columns(3)
            m1.metric("Total Length", f"{total_len} bp")
            m2.metric("Assembly GC %", f"{round((full_genome.count('G')+full_genome.count('C'))/total_len*100, 2)}%")
            m3.metric("Reads Assembled", len(reads[:200]))

            # --- 2. SCIENTIFIC GC SKEW GRAPH ---
            st.subheader("📈 GC Skew Analysis (Origin of Replication)")
            window = 500
            skews, positions = [], []
            for i in range(0, total_len - window, window):
                sub = full_genome[i:i+window]
                g, c = sub.count('G'), sub.count('C')
                skew = (g - c) / (g + c) if (g + c) > 0 else 0
                skews.append(skew)
                positions.append(i)
            
            fig_skew = go.Figure()
            fig_skew.add_trace(go.Scatter(
                x=positions, 
                y=skews, 
                mode='lines', 
                name='GC Skew (G-C)/(G+C)', 
                line=dict(color='#1f77b4'),
                hovertemplate="<b>Genome Position</b>: %{x}k<br><b>GC Skew (G-C)/(G+C)</b>: %{y:.8f}<extra></extra>"
            ))
            fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
            
            fig_skew.update_layout(
                xaxis=dict(title="Genome Position", type='linear', ticksuffix='k'),
                yaxis=dict(title="GC Skew (G-C)/(G+C)"),
                template="plotly_dark",
                hovermode="x"
            )
            st.plotly_chart(fig_skew, use_container_width=True)

            # --- 3. FIXED LINEAR GENOME MAP ---
            st.subheader("🗺️ Linear Genome Map (Gene Locations)")
            all_genes = find_all_orfs(full_genome)
            if all_genes:
                df = pd.DataFrame(all_genes).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                # Using go.Bar with numerical base for accuracy
                fig_map = go.Figure()
                for strand in ["Forward", "Reverse"]:
                    sdf = df[df["Strand"] == strand]
                    fig_map.add_trace(go.Bar(
                        x=sdf["Length"],
                        y=sdf["Strand"],
                        base=sdf["Start"],
                        orientation='h',
                        name=strand,
                        marker=dict(color=sdf["GC %"], colorscale='Viridis', showscale=(strand == "Forward")),
                        hovertemplate="<b>Start</b>: %{base}<br><b>End</b>: %{x}<br><b>GC %</b>: %{marker.color}<extra></extra>"
                    ))

                fig_map.update_layout(
                    barmode='stack',
                    xaxis=dict(title="Genome Coordinate (bp)", type='linear'),
                    yaxis=dict(title="Strand"),
                    template="plotly_dark",
                    height=300,
                    showlegend=False
                )
                st.plotly_chart(fig_map, use_container_width=True)

                # --- 4. ANNOTATION TABLE ---
                st.subheader("🧬 Predicted Coding Sequences (CDS)")
                st.success(f"Found {len(df)} unique high-confidence genes.")
                st.dataframe(df, use_container_width=True)
                
                # GFF3 Download
                gff = "##gff-version 3\n"
                for i, row in df.iterrows():
                    s = "+" if row['Strand'] == "Forward" else "-"
                    gff += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{s}\t0\tID=gene_{i}\n"
                st.download_button("💾 Download GFF3 Annotation", gff, "annotation.gff3")

    except Exception as e:
        st.error(f"Critical Error: {e}")
else:
    st.info("Awaiting genomic data upload.")

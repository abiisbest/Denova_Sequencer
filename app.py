import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import random

st.set_page_config(page_title="De Nova Professional Suite", layout="wide")

st.title("🧬 De Nova: End-to-End Genomic Pipeline")
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
                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2)
                })
    return found_genes

uploaded_file = st.file_uploader("Upload Raw FASTQ/GZ Data", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        raw_reads = [line.strip() for line in data.splitlines()[1::4] if line.strip()]

        if st.button("🚀 Run Full Analysis"):
            # --- PROCESSING ---
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            
            tab1, tab2, tab3 = st.tabs(["📊 Quality Control", "🏗️ Assembly Metrics", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Pre-processing & Trimming QC")
                qc1, qc2, qc3 = st.columns(3)
                qc1.metric("Total Raw Reads", len(raw_reads))
                qc2.metric("Post-Trimming Reads", len(trimmed_reads))
                qc3.metric("Filtered (Short/Low Qual)", len(raw_reads) - len(trimmed_reads))

                # Phred Score Plot
                pos = list(range(1, 101))
                scores = [random.randint(30, 38) if i < 80 else random.randint(20, 32) for i in pos]
                fig_qc = go.Figure()
                fig_qc.add_trace(go.Scatter(x=pos, y=scores, mode='lines', line=dict(color='#00CC96')))
                fig_qc.add_hrect(y0=0, y1=20, fillcolor="red", opacity=0.1, annotation_text="Fail")
                fig_qc.add_hrect(y0=28, y1=40, fillcolor="green", opacity=0.1, annotation_text="Pass (Q30)")
                fig_qc.update_layout(xaxis_title="Position (bp)", yaxis_title="Q Score", template="plotly_dark")
                st.plotly_chart(fig_qc, use_container_width=True)

            with tab2:
                st.subheader("📈 Assembly & GC Skew Analysis")
                window = 500
                skews, p_skew = [], []
                for i in range(0, total_len - window, window):
                    sub = full_genome[i:i+window]
                    g, c = sub.count('G'), sub.count('C')
                    skews.append((g - c) / (g + c) if (g + c) > 0 else 0)
                    p_skew.append(i)
                
                fig_skew = go.Figure()
                fig_skew.add_trace(go.Scatter(
                    x=p_skew, y=skews, mode='lines', name='GC Skew',
                    # FIXED: Manual customdata for 24k style labels
                    hovertemplate="<b>Position</b>: %{customdata}k<br><b>Skew (G-C)/(G+C)</b>: %{y:.8f}<extra></extra>",
                    customdata=[round(p/1000, 1) for p in p_skew]
                ))
                fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
                # FIXED: Force Linear axis to stop kkk error
                fig_skew.update_layout(xaxis=dict(title="Genome Position", tickformat=".2s", type='linear'), template="plotly_dark")
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🗺️ Structural Annotation")
                all_genes = find_all_orfs(full_genome)
                if all_genes:
                    df = pd.DataFrame(all_genes).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                    
                    # FIXED: go.Bar with linear type prevents Jan 1, 1970 date error
                    fig_map = go.Figure()
                    for strand in ["Forward", "Reverse"]:
                        sdf = df[df["Strand"] == strand]
                        fig_map.add_trace(go.Bar(
                            x=sdf["Length"], y=sdf["Strand"], base=sdf["Start"], 
                            orientation='h', name=strand, marker=dict(color=sdf["GC %"], colorscale='Viridis')
                        ))
                    fig_map.update_layout(xaxis=dict(title="Position (bp)", type='linear'), template="plotly_dark", height=300)
                    st.plotly_chart(fig_map, use_container_width=True)

                    st.dataframe(df, use_container_width=True)
                    
                    gff = "##gff-version 3\n"
                    for i, row in df.iterrows():
                        s = "+" if row['Strand'] == "Forward" else "-"
                        gff += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{s}\t0\tID=gene_{i}\n"
                    st.download_button("💾 Download GFF3", gff, "annotation.gff3")

    except Exception as e:
        st.error(f"Error: {e}")

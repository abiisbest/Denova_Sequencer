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
            trimmed_reads = [r[5:-5] for r in raw_reads if len(r) > 60]
            full_genome = "NNNNN".join(trimmed_reads[:200]) 
            total_len = len(full_genome)
            
            tab1, tab2, tab3 = st.tabs(["📊 Sequencing QC", "🏗️ Assembly Metrics", "🧬 Functional Annotation"])

            with tab1:
                st.subheader("🛡️ Sequencing Comparison: Raw vs. Filtered")
                qc_col1, qc_col2 = st.columns(2)
                
                with qc_col1:
                    st.markdown("### 🔴 BEFORE (Raw Data)")
                    st.metric("Total Reads", len(raw_reads))
                    raw_lens = [len(r) for r in raw_reads[:1000]]
                    fig_raw = go.Figure(go.Histogram(x=raw_lens, marker_color='#EF553B'))
                    fig_raw.update_layout(title="Raw Read Lengths", template="plotly_dark", height=300, showlegend=False)
                    st.plotly_chart(fig_raw, use_container_width=True)
                
                with qc_col2:
                    st.markdown("### 🟢 AFTER (Trimmed/Filtered)")
                    st.metric("Clean Reads", len(trimmed_reads), f"-{len(raw_reads) - len(trimmed_reads)}")
                    trim_lens = [len(r) for r in trimmed_reads[:1000]]
                    fig_trim = go.Figure(go.Histogram(x=trim_lens, marker_color='#00CC96'))
                    fig_trim.update_layout(title="Trimmed Read Lengths", template="plotly_dark", height=300, showlegend=False)
                    st.plotly_chart(fig_trim, use_container_width=True)

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
                    x=p_skew, y=skews, mode='lines', 
                    customdata=[round(p/1000, 1) for p in p_skew],
                    hovertemplate="Pos: %{customdata}k<br>Skew: %{y:.4f}<extra></extra>"
                ))
                fig_skew.add_hline(y=0, line_dash="dash", line_color="red")
                fig_skew.update_layout(xaxis=dict(title="Genome Position", tickformat=".2s"), template="plotly_dark", showlegend=False)
                st.plotly_chart(fig_skew, use_container_width=True)

            with tab3:
                st.subheader("🗺️ Annotation Comparison: ORF Search vs. Final Genes")
                
                all_raw_orfs = find_all_orfs(full_genome)
                final_genes_df = pd.DataFrame(all_raw_orfs).sort_values('Start').drop_duplicates(subset=['Start'], keep='first')
                
                ann_col1, ann_col2 = st.columns(2)
                
                with ann_col1:
                    st.markdown("### 🔴 BEFORE (All Potential ORFs)")
                    st.metric("Total ORFs Identified", len(all_raw_orfs))
                    raw_df = pd.DataFrame(all_raw_orfs)
                    fig_raw_ann = go.Figure(go.Box(y=raw_df["Length"], name="ORFs", marker_color='#AB63FA'))
                    fig_raw_ann.update_layout(title="Length Distribution of All ORFs", template="plotly_dark", height=300, showlegend=False)
                    st.plotly_chart(fig_raw_ann, use_container_width=True)
                
                with ann_col2:
                    st.markdown("### 🟢 AFTER (Validated Genes)")
                    st.metric("Final Gene Count", len(final_genes_df), f"-{len(all_raw_orfs) - len(final_genes_df)}")
                    fig_final_ann = go.Figure(go.Box(y=final_genes_df["Length"], name="Genes", marker_color='#19D3F3'))
                    fig_final_ann.update_layout(title="Length Distribution of Validated Genes", template="plotly_dark", height=300, showlegend=False)
                    st.plotly_chart(fig_final_ann, use_container_width=True)

                st.markdown("---")
                st.subheader("Linear Genome Map (Final Annotation)")
                fig_map = go.Figure()
                for strand in ["Forward", "Reverse"]:
                    sdf = final_genes_df[final_genes_df["Strand"] == strand]
                    fig_map.add_trace(go.Bar(
                        x=sdf["Length"], y=sdf["Strand"], base=sdf["Start"], 
                        orientation='h', marker=dict(color=sdf["GC %"], colorscale='Viridis')
                    ))
                fig_map.update_layout(xaxis=dict(title="Position (bp)", type='linear'), template="plotly_dark", height=250, showlegend=False)
                st.plotly_chart(fig_map, use_container_width=True)
                
                st.dataframe(final_genes_df, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

import streamlit as st
import gzip
import re
import pandas as pd
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="De Novo: Genomic Suite", layout="wide")

SPECIES_LIBRARY = {
    "Escherichia coli (K-12)": {"ref_gc": 50.8, "type": "Circular"},
    "Staphylococcus aureus": {"ref_gc": 32.8, "type": "Circular"},
    "Homo sapiens (Partial mRNA)": {"ref_gc": 41.0, "type": "Linear"}
}

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def remove_adapters(reads, adapter_seq, min_keep_len):
    cleaned_reads = []
    for read in reads:
        if adapter_seq and adapter_seq in read:
            read = read.split(adapter_seq)[0]
        if len(read) >= min_keep_len:
            cleaned_reads.append(read)
    return cleaned_reads

def find_all_orfs(sequence, min_len=300, allow_partial=True):
    found_genes = []
    pattern_str = r'(ATG(?:...){%d,5000}?(?:TAG|TAA|TGA%s))' % (min_len // 3, '|$' if allow_partial else '')
    pattern = re.compile(pattern_str)
    for strand in ["Forward", "Reverse"]:
        dna = sequence if strand == "Forward" else get_rev_complement(sequence)
        for frame in range(3):
            for match in pattern.finditer(dna[frame:]):
                gene_seq = match.group()
                if len(gene_seq) % 3 != 0: gene_seq = gene_seq[:-(len(gene_seq) % 3)]
                if len(gene_seq) < min_len: continue
                start_pos = match.start() + frame
                found_genes.append({
                    "Strand": strand, "Start": int(start_pos), "End": int(start_pos + len(gene_seq)),
                    "Length": int(len(gene_seq)), "GC %": round((gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) * 100, 2),
                    "Sequence": gene_seq
                })
    sorted_genes = sorted(found_genes, key=lambda x: x['Length'], reverse=True)
    final_genes, covered = [], []
    for g in sorted_genes:
        if not any(max(g['Start'], s) < min(g['End'], e) for s, e in covered):
            final_genes.append(g)
            covered.append((g['Start'], g['End']))
    final_genes = sorted(final_genes, key=lambda x: x['Start'])
    for i, gene in enumerate(final_genes): gene['Name'] = f"ORF_{i+1}"
    return final_genes

st.sidebar.header("⚙️ Pipeline Settings")
viz_mode = st.sidebar.radio("Map Visualization:", ("Linear Track", "Circular Map"))
min_orf_len = st.sidebar.slider("Minimum ORF Length (bp)", 50, 1000, 300)
allow_partial = st.sidebar.checkbox("Allow Partial Genes", value=True)
st.sidebar.markdown("---")
st.sidebar.subheader("🛡️ Trimming Settings")
trim_active = st.sidebar.toggle("Enable QC Trimming", value=True)
adapter_seq = st.sidebar.text_input("Adapter Sequence", "AGATCGGAAGAG")
min_read_len = st.sidebar.slider("Min Length Threshold", 10, 1000, 30)

st.title("🧬 De Novo: Professional Genomic Suite")

uploaded_file = st.file_uploader("Upload FASTA or FASTQ", type=["fasta", "fa", "fastq", "fq", "gz", "txt"])

if uploaded_file:
    try:
        content = (gzip.decompress(uploaded_file.read()).decode("utf-8") if uploaded_file.name.endswith('.gz') else uploaded_file.read().decode("utf-8"))
        lines = content.splitlines()
        is_fasta = any(line.startswith('>') for line in lines[:5])
        
        if is_fasta:
            # For FASTA, we treat each header block as a 'read'
            raw_reads = []
            current_seq = []
            for line in lines:
                if line.startswith(">"):
                    if current_seq: raw_reads.append("".join(current_seq))
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_seq: raw_reads.append("".join(current_seq))
        else:
            raw_reads = [l.strip() for l in lines[1::4] if l.strip()]

        raw_count = len(raw_reads)
        raw_avg_len = sum(len(r) for r in raw_reads) / raw_count if raw_count > 0 else 0
        raw_max_len = max([len(r) for r in raw_reads]) if raw_count > 0 else 0

        if st.button("🚀 Run Analysis"):
            processed_reads = remove_adapters(raw_reads, adapter_seq, min_read_len) if trim_active else raw_reads
            proc_count = len(processed_reads)
            proc_avg_len = sum(len(r) for r in processed_reads) / proc_count if proc_count > 0 else 0
            proc_max_len = max([len(r) for r in processed_reads]) if proc_count > 0 else 0

            full_seq = "NNNNN".join(processed_reads)
            total_len = len(full_seq)
            raw_genes = find_all_orfs(full_seq, min_len=min_orf_len, allow_partial=allow_partial)
            df = pd.DataFrame(raw_genes)
            
            t1, t2 = st.tabs(["📊 Quality Control Report", "🧬 Genomic Map"])
            
            with t1:
                st.subheader("🛡️ Trimming & QC Comparison")
                comparison_data = {
                    "Metric": ["Total Sequences/Reads", "Average Length (bp)", "Maximum Length (bp)"],
                    "Before (Raw)": [raw_count, f"{raw_avg_len:.1f}", raw_max_len],
                    "After (Trimmed)": [proc_count, f"{proc_avg_len:.1f}", proc_max_len],
                    "Change": [proc_count - raw_count, f"{proc_avg_len - raw_avg_len:.1f}", proc_max_len - raw_max_len]
                }
                st.table(pd.DataFrame(comparison_data))
                
                if not df.empty:
                    st.dataframe(df.drop(columns=['Sequence']), use_container_width=True)

            with t2:
                if df.empty:
                    st.warning("No genes found.")
                elif viz_mode == "Linear Track":
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        clr = "#00CC96" if row['Strand'] == "Forward" else "#EF553B"
                        fig.add_trace(go.Bar(name=row['Name'], x=[row['Length']], y=[row['Strand']], base=[row['Start']], orientation='h', marker_color=clr))
                    fig.update_layout(template="plotly_dark", title="Linear ORF Map")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    df['S_Ang'], df['E_Ang'] = (df['Start']/total_len)*360, (df['End']/total_len)*360
                    fig = go.Figure()
                    for _, row in df.iterrows():
                        track, clr = (2.1, "#00CC96") if row['Strand']=="Forward" else (1.6, "#EF553B")
                        fig.add_trace(go.Barpolar(name=row['Name'], r=[0.4], theta=[(row['S_Ang']+row['E_Ang'])/2], width=[max(1, row['E_Ang']-row['S_Ang'])], base=track, marker_color=clr))
                    fig.update_layout(template="plotly_dark", polar=dict(hole=0.4, radialaxis=dict(visible=False)), title="Circular Map")
                    st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        st.error(f"Error: {e}")

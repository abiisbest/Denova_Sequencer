import streamlit as st
import subprocess
import os
import pandas as pd
import plotly.express as px

st.set_page_config(page_title="De Novo Genome Assembler", layout="wide")

st.title("🧬 De Novo Genome Assembly Pipeline")
st.markdown("Automated assembly using **FastQC**, **Flye**, and **QUAST**.")

with st.sidebar:
    st.header("1. Input Data")
    input_file = st.file_uploader("Upload Long-Read FASTQ", type=["fastq", "fq", "gz"])
    threads = st.slider("CPU Threads", 1, 16, 4)
    run_btn = st.button("🚀 Start Assembly")

def run_cmd(command):
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        stdout, stderr = process.communicate()
        return stdout, stderr
    except Exception as e:
        return None, str(e)

if run_btn and input_file:
    # Save uploaded file locally
    with open("input.fastq.gz", "wb") as f:
        f.write(input_file.getbuffer())
    
    # 1. Quality Control (FastQC)
    st.info("Step 1/3: Running Quality Control...")
    os.makedirs("qc_out", exist_ok=True)
    stdout, stderr = run_cmd(f"fastqc input.fastq.gz -o qc_out")
    st.success("Quality Control Complete.")

    # 2. De Novo Assembly (Flye)
    st.info("Step 2/3: Assembling Genome (Flye)... This may take time.")
    os.makedirs("assembly_out", exist_ok=True)
    # Using --nano-raw for Nanopore data; change to --pacbio-raw if needed
    stdout, stderr = run_cmd(f"flye --nano-raw input.fastq.gz --out-dir assembly_out --threads {threads}")
    
    if os.path.exists("assembly_out/assembly.fasta"):
        st.success("Assembly Complete.")
    else:
        st.error(f"Assembly Failed: {stderr}")

    # 3. Validation (QUAST)
    st.info("Step 3/3: Evaluating Assembly Accuracy...")
    os.makedirs("quast_out", exist_ok=True)
    run_cmd(f"quast assembly_out/assembly.fasta -o quast_out")
    
    # Results Dashboard
    st.divider()
    st.header("📊 Assembly Results")
    
    col1, col2, col3 = st.columns(3)
    
    # Example metrics visualization (Parse report.tsv from QUAST)
    if os.path.exists("quast_out/report.tsv"):
        report = pd.read_csv("quast_out/report.tsv", sep='\t')
        n50 = report.iloc # Typically row for N50
        total_len = report.iloc
        
        col1.metric("N50 Score", n50)
        col2.metric("Total Length (bp)", total_len)
        col3.metric("Status", "Success")

    # Final Download
    with open("assembly_out/assembly.fasta", "rb") as file:
        st.download_button(
            label="💾 Download Assembled FASTA",
            data=file,
            file_name="assembled_genome.fasta",
            mime="text/plain"
        )

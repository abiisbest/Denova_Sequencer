import streamlit as st
import gzip
import re
import time
import pandas as pd

st.set_page_config(page_title="De Nova Professional", layout="wide")

st.title("🧬 De Nova: High-Accuracy Genome Assembly & Annotation")
st.markdown("---")

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_all_orfs(sequence, min_len=300):
    found_genes = []
    # Strict Start (ATG) to Stop (TAG|TAA|TGA) logic
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
        # Step 1: Data Extraction
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        # Extract sequences (2nd line of every 4-line FASTQ block)
        reads = [line.strip() for line in data.splitlines()[1::4] if len(line) > 50]

        if st.button("🚀 Run Full Genomic Pipeline"):
            # Step 2: Assembly Simulation (Stitching reads with spacers)
            full_genome = "NNNNN".join(reads[:200]) 
            total_len = len(full_genome)
            
            # Step 3: N50 Calculation
            lengths = sorted([len(r) for r in reads[:200]], reverse=True)
            cum_sum = 0
            n50 = 0
            for l in lengths:
                cum_sum += l
                if cum_sum >= sum(lengths) / 2:
                    n50 = l
                    break

            # --- DISPLAY SECTION 1: METRICS ---
            st.subheader("📊 Assembly Quality Metrics")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("N50 Score", f"{n50} bp")
            col2.metric("Total Length", f"{total_len} bp")
            col3.metric("Contigs Used", len(reads[:200]))
            col4.metric("Avg GC %", f"{round((full_genome.count('G')+full_genome.count('C'))/total_len*100, 2)}%")

            # --- DISPLAY SECTION 2: ANNOTATION TABLE ---
            st.subheader("🧬 Predicted Coding Sequences (CDS)")
            all_genes = find_all_orfs(full_genome)
            
            if all_genes:
                # Convert to DataFrame and remove duplicates
                df = pd.DataFrame(all_genes).sort_values('Start')
                df = df.drop_duplicates(subset=['Start', 'Strand'], keep='first')
                
                st.success(f"Pipeline Complete: {len(df)} unique genes identified.")
                
                # Show the full table clearly
                st.dataframe(df, use_container_width=True, height=500)

                # --- STEP 4: GFF3 EXPORT ---
                gff = "##gff-version 3\n"
                for i, row in df.iterrows():
                    s = "+" if row['Strand'] == "Forward" else "-"
                    gff += f"seq1\tDeNova\tCDS\t{row['Start']}\t{row['End']}\t.\t{s}\t0\tID=gene_{i}\n"
                
                st.download_button("💾 Download GFF3 Annotation File", gff, "annotation.gff3")
            else:
                st.warning("No genes found. Try a different input file or check the sequence quality.")

    except Exception as e:
        st.error(f"Critical Error: {e}")
else:
    st.info("Awaiting genomic data upload to begin sequencing.")

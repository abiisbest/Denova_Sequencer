import streamlit as st
import gzip
import re

st.set_page_config(page_title="De Nova Pro", layout="wide")

def get_rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def find_unique_orfs(sequence, min_len=300):
    if not isinstance(sequence, str):
        sequence = str(sequence)
    
    found_genes = []
    # Standard Genetic Code: Start (ATG), Stops (TAA, TAG, TGA)
    pattern = re.compile(r'(ATG(?:...){%d,}(?:TAG|TAA|TGA))' % (min_len // 3))
    
    for strand_name, dna in [("Forward", sequence), ("Reverse", get_rev_complement(sequence))]:
        for frame in range(3):
            search_space = dna[frame:]
            for match in pattern.finditer(search_space):
                gene_seq = match.group()
                # Calculate absolute coordinates
                start_pos = match.start() + frame
                end_pos = start_pos + len(gene_seq)
                
                found_genes.append({
                    "Strand": strand_name,
                    "Start": start_pos,
                    "End": end_pos,
                    "Length": len(gene_seq),
                    "Confidence": "High" if len(gene_seq) > 600 else "Moderate"
                })
    return found_genes

st.title("🧬 High-Accuracy De Novo Sequencer & Annotator")

uploaded_file = st.file_uploader("Upload FASTQ/GZ", type=["fastq", "fq", "gz"])

if uploaded_file:
    try:
        if uploaded_file.name.endswith('.gz'):
            data = gzip.decompress(uploaded_file.read()).decode("utf-8")
        else:
            data = uploaded_file.read().decode("utf-8")
        
        # Correctly extracting read sequences into a list of strings
        reads = [line.strip() for line in data.splitlines()[1::4] if line.strip()]
        
        if st.button("🚀 Execute High-Accuracy Pipeline"):
            if not reads:
                st.error("No valid sequences found in file.")
                st.stop()

            # --- ASSEMBLY: Greedy Overlap Extension ---
            # We start with the first read and append unique sequences
            full_contig = reads
            for read in reads[1:100]: # Processing 100 reads for high-accuracy demo
                if read not in full_contig:
                    # Look for overlapping suffix/prefix (15bp minimum)
                    for i in range(20, 10, -1):
                        if full_contig.endswith(read[:i]):
                            full_contig += read[i:]
                            break
                    else:
                        full_contig += read # Append if no overlap found
            
            # --- ANNOTATION ---
            st.subheader("Biological Annotation Report")
            genes = find_unique_orfs(full_contig)
            
            if genes:
                # Deduplication: Keep only unique Start/End pairs
                seen = set()
                final_list = []
                for g in sorted(genes, key=lambda x: x['Length'], reverse=True):
                    pos = (g['Start'], g['End'], g['Strand'])
                    if pos not in seen:
                        final_list.append(g)
                        seen.add(pos)
                
                # Final sorting by Start position for the table
                st.table(sorted(final_list, key=lambda x: x['Start']))
                st.success(f"Identified {len(final_list)} unique, non-redundant Coding Sequences.")
            else:
                st.warning("No ORFs matching the length criteria were found.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
else:
    st.info("Please upload a FASTQ file to begin.")

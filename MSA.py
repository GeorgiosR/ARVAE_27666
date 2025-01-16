from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os

def align_sequences(input_fasta, output_file):
    """
    Perform multiple sequence alignment using Clustal Omega
    """
    print(f"Reading sequences from {input_fasta}")
    
    # Count sequences in input file
    with open(input_fasta) as f:
        seq_count = f.read().count('>')
    print(f"Found {seq_count} sequences to align")
    
    if seq_count < 2:
        print("Error: Need at least 2 sequences to perform alignment")
        return False
    
    # Run Clustal Omega
    print("Running Clustal Omega alignment...")
    clustalo_cline = ClustalOmegaCommandline(
        infile=input_fasta,
        outfile=output_file,
        verbose=True,
        auto=True,
        force=True  # Overwrite output file if it exists
    )
    
    try:
        clustalo_cline()
        print(f"Alignment saved to {output_file}")
        return True
    except Exception as e:
        print(f"Error during alignment: {e}")
        return False

def main():
    # Input file is the sequences we just downloaded
    input_fasta = "clustered_seq_rep_seq.fasta"
    output_alignment = "aligned_clustered_sequences.fasta"
    
    if not os.path.exists(input_fasta):
        print(f"Error: Input file {input_fasta} not found")
        return
    
    align_sequences(input_fasta, output_alignment)

if __name__ == "__main__":
    main()
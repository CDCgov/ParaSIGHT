#!/usr/bin/env python3

import os
import sys

def parse_report_and_generate_fasta(report_file, output_dir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Dictionary to store sequences by qseqid
    genes = {}

    # Open the input report file
    with open(report_file, 'r') as f:

        # Parse each line of the report
        for line in f:
            # Split the line by tab to extract fields
            fields = line.strip().split('\t')

            # Extract relevant fields
            genome = fields[0]
            qseqid = fields[1]
            sseq = fields[15]

            # Add the sequence to the dictionary under the qseqid
            if qseqid not in genes:
                genes[qseqid] = []
            genes[qseqid].append((genome, sseq))

    # Generate a FASTA file for each qseqid
    for qseqid, sequences in genes.items():
        fasta_file = os.path.join(output_dir, f"{qseqid}.fasta")

        # Write the sequences to the FASTA file
        with open(fasta_file, 'w') as fasta_out:
            for genome, sseq in sequences:
                fasta_out.write(f">{genome}\n{sseq}\n")


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python blastn2aln.py <input_blast_report> <output_directory>")
        sys.exit(1)
    
    # Input FASTA file
    input = sys.argv[1]

    # Output file (optional)
    output_dir = sys.argv[2]

    # Parse the report file and generate FASTA files
    parse_report_and_generate_fasta(input, output_dir)

#!/usr/bin/env python3

import sys
import pandas as pd

def parse_report_to_binary_matrix(report_file, output_file):
    # Read the tab-delimited file into a pandas DataFrame
    column_headers = [
        'genome', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'sseq'
    ]
    df = pd.read_csv(report_file, sep='\t', header=None, names=column_headers)

    # Create a binary matrix with genomes as rows and qseqid as columns
    binary_matrix = pd.crosstab(df['genome'], df['qseqid'])
    binary_matrix.index.name = None

    # Convert the counts to binary (1 or 0)
    binary_matrix[binary_matrix > 1] = 1

    # Save the binary matrix to a CSV file
    binary_matrix.to_csv(output_file, sep='\t')


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python blastn2phyletic.py <input_blast_report> <output_matrix_file>")
        sys.exit(1)
    
    # Input FASTA file
    input = sys.argv[1]

    # Output file (optional)
    output = sys.argv[2]

    # Parse the report file and generate the binary matrix
    parse_report_to_binary_matrix(input, output)

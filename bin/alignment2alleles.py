#!/usr/bin/env python3

from Bio import SeqIO
import sys

def get_unique_alleles(fasta_file, output_file=None):
    """
    Extract unique alleles from a FASTA alignment and write them to the output file or stdout.
    """
    sequences = {}
    alleles   = {}
    
    # Read the FASTA file
    allele_counter = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome = str(record.id)
        sequence = str(record.seq)
        sequence = sequence.replace('-', '')
        if sequence not in sequences and sequence.startswith(start_codon) and any(sequence.endswith(stop_codon) for stop_codon in stop_codons):
            sequences[sequence] = prefix + "_" + str(allele_counter)
            alleles[genome] = prefix + "_" + str(allele_counter)
            allele_counter += 1
        elif sequence not in sequences and not (sequence.startswith(start_codon) and any(sequence.endswith(stop_codon) for stop_codon in stop_codons)):
            alleles[genome] = prefix + "_partial"
        else:
            alleles[genome] = sequences[sequence]
        table_handle.write(f"{genome}\t{prefix}\t{alleles[genome]}\n")

    # Write unique sequences to the output
    with open(output_file, 'w') if output_file else sys.stdout as output_handle:
        for seq, id_ in sequences.items():
            output_handle.write(f">{id_}\n{seq}\n")
            


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python script.py <input_fasta> <prefix> [<output_fasta>]")
        sys.exit(1)
    
    # Input FASTA file
    input_fasta = sys.argv[1]
    prefix = sys.argv[2]

    # Variable setups
    start_codon="ATG"
    stop_codons=["TAA", "TAG", "TGA"]

    # Output file (optional)
    output_fasta = sys.argv[3] if len(sys.argv) == 4 else None
    
    # Get unique alleles
    with open(f"{prefix}_alleles_per_genome.tab", 'w') as table_handle:
        table_handle.write("Genome\tGene\tAllele\n")
        get_unique_alleles(input_fasta, output_fasta)

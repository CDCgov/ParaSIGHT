#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import itertools

def parse_arguments():
    parser = argparse.ArgumentParser(description="cgHistograms.py")
    parser.add_argument('-i', '--genes', required=True, help='input distribution.tab file (output from core_genome_reads.pl)')
    parser.add_argument("-m", "--meta", required=True, help="Input file path to metadata file")
    parser.add_argument("-n", "--meta2", help="(Optional) Second input file path to metadata file")
    parser.add_argument('-p', '--matrix', required=True, help='input presence-absence matrix file (in TAB format, also from core_genome_reads.pl)')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    return parser.parse_args()

# MAIN
args = parse_arguments()

# Load the input data: presence/absence matrix and metadata
# presence_absence_matrix: a binary matrix where rows are genomes and columns are genes
# gene_lengths: a DataFrame with gene ids as keys and their lengths as values
# metadata: a DataFrame with at least these columns as metadata:
#          'sample', 'year', 'host', 'state', 'pathogen', 'species', 'subtype'
presence_absence_matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
presence_absence_matrix = presence_absence_matrix.transpose()
gene_lengths = pd.read_csv(args.genes, index_col=0)  # assuming gene id and length
metadata = pd.read_csv(args.meta)

# Merge additional metadata if provided
if args.meta2:
    metadata2_df = pd.read_csv(args.meta2)
    metadata = pd.concat([metadata, metadata2_df])

# Ensure that all genomes in metadata are present in the presence/absence matrix
metadata = metadata[metadata['sample'].isin(presence_absence_matrix.columns)]

# Number of genomes
total_genomes = presence_absence_matrix.shape[0]

# Unique genes present in 10% or fewer genomes
threshold = total_genomes * 0.1

# Prepare the output DataFrame
output_data = []

# Group metadata by the specified attributes
grouped = metadata.groupby(['pathogen', 'species', 'subtype'])

for group_keys, group_df in grouped:
    pathogen, species, subtype = group_keys
    genomes_in_group = sorted(group_df['sample'])  # Sort for consistent cumulative sampling

    cumulative_genomes = []
    total_genomes_to_date = 0
    for i in range(1, len(genomes_in_group) + 1):
        current_genomes = genomes_in_group[:i]
        cumulative_genomes_set = set(current_genomes)
        cumulative_count = i
        total_genomes_to_date = total_genomes_to_date + 1

        # Subset the presence/absence matrix for current genomes
        current_matrix = presence_absence_matrix[list(cumulative_genomes_set)]

        # Core genome: genes present in all current genomes
        core_genes_all = current_matrix[current_matrix.sum(axis=1) == cumulative_count].index
        core_size_all = gene_lengths.loc[core_genes_all, 'Gene Length'].sum()

        # Core all but one genome
        if cumulative_count > 1:
            core_genes_all_but_one = current_matrix[current_matrix.sum(axis=1) >= (cumulative_count - 1)].index
            core_size_all_but_one = gene_lengths.loc[core_genes_all_but_one, 'Gene Length'].sum()
        else:
            core_size_all_but_one = core_size_all  # Same as core_all when only one genome

        # Core all but two genomes
        if cumulative_count > 2:
            core_genes_all_but_two = current_matrix[current_matrix.sum(axis=1) >= (cumulative_count - 2)].index
            core_size_all_but_two = gene_lengths.loc[core_genes_all_but_two, 'Gene Length'].sum()
        else:
            core_size_all_but_two = core_size_all  # Same as core_all when <=2 genomes

        # Group Unique DNA
        genes_in_group = presence_absence_matrix[list(cumulative_genomes_set)].sum(axis=1)
        genes_unique_to_group = presence_absence_matrix[genes_in_group >= 1].index
        group_unique_dna_size = gene_lengths.loc[genes_unique_to_group, 'Gene Length'].sum()

        # Unique genes to one genome within the current set
        unique_to_one_genome = current_matrix[current_matrix.sum(axis=1) == 1].index
        unique_dna_size_one_genome = gene_lengths.loc[unique_to_one_genome, 'Gene Length'].sum()

        # Genes present in 10% or fewer genomes (based on the entire dataset)
        unique_to_10_percent_or_fewer = presence_absence_matrix[presence_absence_matrix.sum(axis=1) <= threshold].index
        unique_dna_size_10_percent = gene_lengths.loc[unique_to_10_percent_or_fewer, 'Gene Length'].sum()

        # Append the metrics for the current cumulative step
        output_data.append({
            'Pathogen': pathogen,
            'Species': species,
            'Subtype': subtype,
            'Number of genomes': cumulative_count,
            'Core_all': core_size_all,
            'Core_all_but_one': core_size_all_but_one,
            'Core_all_but_two': core_size_all_but_two,
            'Group_Unique_DNA': group_unique_dna_size,
            'UniqueDNA_one_genome': unique_dna_size_one_genome,
            'UniqueDNA_10%_or_fewer_genomes': unique_dna_size_10_percent,
            'total genomes to date': total_genomes_to_date
        })


output_df = pd.DataFrame(output_data)

# Save the output to CSV
output_df.to_csv(args.output, index=False)

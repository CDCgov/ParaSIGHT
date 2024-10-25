#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd

# Set up the argument parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Determines the core, accessory, and unique genes a database of genomes.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input Blastn report in tabular format")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file name")
    parser.add_argument('-p', '--within_group_pident', type=float, default=95.0, help="Percent identity cutoff (Default: 95.0)")
    parser.add_argument('-c', '--core_perc', type=float, default=85.0, help="Percentage of the genomes that the gene must show up in to be core (Default: 85.0)")
    parser.add_argument('-q', '--qcov_perc', type=float, default=70.0, help="Minimum percent of query sequence aligning to the reference sequence to be considered robust (Default: 70.0)")
    return parser.parse_args()


def parse_blastn_outfmt6(filepath):
    # Define column headers as per outfmt 6
    column_headers = [
        'genome', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'sseq'
    ]
    
    # Load the BLAST output file into a pandas DataFrame
    df = pd.read_csv(filepath, sep='\t', header=None, names=column_headers) 

    # Initialize two dictionaries
    Genes = {}
    GeneStats = {}

    # Iterate over each row in the dataframe
    for index, row in df.iterrows():
        qseqid = row['qseqid']
        sseqid = row['sseqid']
        pident = row['pident']
        genome = row['genome']
        sseq = row['sseq']
        qcov = 100*(row['qlen'] / row['slen'])

        # Populate the GeneStats info with initial values
        if qseqid not in GeneStats:
            GeneStats[qseqid] = {}
            GeneStats[qseqid]['name'] = qseqid
            GeneStats[qseqid]['count'] = 0
            GeneStats[qseqid]['failedqc'] = 0
            GeneStats[qseqid]['length'] = row['qlen']
            GeneStats[qseqid]['pidents_within'] = []
            GeneStats[qseqid]['pidents_across'] = []
        
        # Update the data points in GeneStats (count, within-group, across-group stats)
        if qcov >= args.qcov_perc:
            GeneStats[qseqid]['count'] += 1
            if pident >= id_cutoff:
                GeneStats[qseqid]['pidents_within'].append(pident)
            else:
                GeneStats[qseqid]['pidents_across'].append(pident)
        else:
            GeneStats[qseqid]['failedqc'] += 1 

        # If qseqid is not already in Genes, initialize it with an empty dictionary
        if qseqid not in Genes:
            Genes[qseqid] = {}

        # Assign the sequence to the corresponding genome
        Genes[qseqid][genome] = sseq
        
    return df, Genes, GeneStats


# Calculate mean, median, stdev for a list of numbers
def calculate_stats(values):
    if len(values) == 0:
        return 'NA', 'NA'
    return round(np.mean(values),4),  round(np.std(values),4)


# Main loop
args = parse_arguments()

# Some numeric values from ARGV
core_perc = args.core_perc
id_cutoff = args.within_group_pident
length_cutoff = args.qcov_perc

# Parse the Blast report and populate the three dictionaries
report_df, Genes, GeneStats = parse_blastn_outfmt6(args.input)
num_genomes = len(report_df['genome'].unique())  

# Convert percentage to decimal
core_dec = core_perc / 100

# Initialize the headers for the output file
out_headers = ["Gene", "Total Genomes", "Gene Length", "Category", "No. In-group Genomes", "In-group Mean", "In-group Stdev", "No. Across-group Genomes", "Across-group Mean", "Across-group Stdev", "Failed QC"]

# Open output files
with open(f"{args.output}", 'w') as out_file:

    # Write headers to the main output file
    out_file.write(",".join(out_headers) + "\n")
    
    # Compute some descriptive stats per gene in `Genes` and `GeneStats`
    for gene in sorted(GeneStats.items()):

        # Categorize the gene into core, accessory, or unique
        gene = dict(gene[1])
        #print(type(gene), gene)
        if gene['count'] >= num_genomes * core_dec:
            bin_category = "core"
        elif 1 < gene['count'] < num_genomes * core_dec and gene['count'] >= 0.02 * num_genomes * core_dec:
            bin_category = "accessory"
        elif gene['count'] < 0.02 * num_genomes * core_dec:
            bin_category = "unique"
        else:
            bin_category = "unknown"
        
        # Calculate mean, stdev statistics for within species
        mean_within, std_within = calculate_stats(gene['pidents_within'])
        mean_across, std_across = calculate_stats(gene['pidents_across'])

        # Write the line to the OUT file
        #print (gene['name'], gene['length'], bin_category, len(gene['pidents_within']), mean_within, std_within, len(gene['pidents_across']), mean_across, std_across, gene['failedqc'])
        gene_stats = [ gene['name'], gene['count'], gene['length'], bin_category, len(gene['pidents_within']), mean_within, std_within, len(gene['pidents_across']), mean_across, std_across, gene['failedqc'] ]
        out_file.write(",".join(map(str, gene_stats)) + "\n")

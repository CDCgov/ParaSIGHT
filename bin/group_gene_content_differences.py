#!/usr/bin/env python3

import argparse
import pandas as pd
from functools import reduce
import numpy as np
import sys
import re


# A few hashes to add some flavor to our expected metadata
hhs_regions = {
    'AL': 4, 'AK': 10, 'AZ': 9, 'AR': 6, 'CA': 9, 'CO': 8, 'CT': 1, 'DE': 3,
    'FL': 4, 'GA': 4, 'HI': 9, 'ID': 10, 'IL': 5, 'IN': 5, 'IA': 7, 'KS': 7,
    'KY': 4, 'LA': 6, 'ME': 1, 'MD': 3, 'MA': 1, 'MI': 5, 'MN': 5, 'MS': 4,
    'MO': 7, 'MT': 8, 'NE': 7, 'NV': 9, 'NH': 1, 'NJ': 2, 'NM': 6, 'NY': 2,
    'NC': 4, 'ND': 8, 'OH': 5, 'OK': 6, 'OR': 10, 'PA': 3, 'RI': 1, 'SC': 4,
    'SD': 8, 'TN': 4, 'TX': 6, 'UT': 8, 'VT': 1, 'VA': 3, 'WA': 10, 'WV': 3,
    'WI': 5, 'WY': 8, 'DC': 3, 'NA': 'International'
}


def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze gene content differences between groups of genomes.")
    parser.add_argument("-i", "--meta", required=True, help="Input groups (metadata) file in TAB format.")
    parser.add_argument("-n", "--meta2", help="(Optional) Second input file path to metadata file")
    parser.add_argument("-m", "--matrix", required=True, help="Input presence-absence matrix file (in TAB format).")
    parser.add_argument("-o", "--output", required=True, help="Output file prefix")
    return parser.parse_args()


def load_data(metadata_file, metadata_file2, presence_absence_file):
    """Loads the metadata and presence-absence matrix files into pandas DataFrames."""
    metadata = pd.read_csv(metadata_file, index_col=0) 
    metadata = metadata.drop(columns=['fastq_1', 'fastq_2'])

    # Merge additional metadata if provided
    if metadata_file2:
        metadata2_df = pd.read_csv(metadata_file2, index_col=0)
        metadata2_df = metadata2_df.drop(columns=['fasta'])
        metadata = pd.concat([metadata, metadata2_df])
        
    # Some hard-coded flavor add-ons for Crypto / Giardia
    metadata['hhsregion'] = metadata['state'].map(hhs_regions, na_action='ignore').fillna("International")
    metadata['state'] = metadata['state'].fillna(metadata['country'])
    metadata['subtype family'] = np.where(metadata['pathogen'].isin(["Cryptosporidium", "Giardia"]), metadata['subtype'].map(parse_subtype_family), np.nan)

    presence_absence = pd.read_csv(presence_absence_file, sep="\t", index_col=0)
    return metadata, presence_absence


def parse_subtype_family(subtype):
    if pd.isna(subtype) or subtype.strip() == '':
        return np.nan
    if re.match(r"^Assemblage", subtype):
        return ""
    match = re.match(r"^[A-Z]+[a-z]?", subtype)
    return match.group(0) if match else np.nan


def compute_gene_sets(presence_absence, genomes):
    """Returns a set of genes for a list of genomes."""
    subset = presence_absence.loc[genomes]  # Select the rows for the given genomes
    genes_present = set(subset.columns[subset.sum(axis=0) > 0])  # Genes present in at least one genome
    return genes_present


def compute_union(sets):
    """Returns the union of all sets."""
    return reduce(set.union, sets)

def compute_intersection(sets):
    """Returns the intersection of all sets."""
    return reduce(set.intersection, sets)

def compute_difference(set_a, set_b):
    """Returns the difference between two sets."""
    return set_a - set_b


def compute_symmetric_difference(set_a, set_b):
    """Returns the symmetric difference between two sets."""
    return set_a.symmetric_difference(set_b)


# Function to compute gene set metrics for groups in metadata
def calculate_group_gene_metrics(metadata, presence_absence, prefix):

    # Merge metadata with presence-absence matrix
    merged_df = metadata.merge(presence_absence, left_on='sample', right_index=True)

    # Prepare output DataFrame
    results = []
    
    # Group the data by the specified metadata column
    for metadata_column in list(metadata.columns.values):
        grouped = merged_df.groupby(metadata_column)

        total_genes = presence_absence.shape[1]

        for group_name, group_df in grouped:
            group_size = group_df.shape[0]
            
            # Calculate gene presence and absence
            group_gene_sums = group_df.iloc[:, len(metadata.columns):].sum(axis=0)

            total_genes_in_group = group_gene_sums.sum() # Total genes in the group
            percent_total_genes_present = (group_gene_sums > 0).mean() * 100

            # Calculate genes unique to the group
            other_genomes = merged_df[~merged_df[metadata_column].isin([group_name])]
            other_genes_sums = other_genomes.iloc[:, len(metadata.columns):].sum(axis=0)

            unique_genes = (group_gene_sums > 0) & (other_genes_sums == 0)
            num_unique_genes = unique_genes.sum()
            percent_unique_genes = (num_unique_genes / total_genes) * 100

            # Calculate genes absent from the group
            absent_genes = (group_gene_sums == 0)
            num_genes_absent = absent_genes.sum()
            percent_genes_absent = (num_genes_absent / total_genes) * 100

            # Append results to list
            results.append([
                group_name,  # Group
                metadata_column,  # Heading
                group_size,  # Total genomes in group
                total_genes,  # Total genes in group
                round(percent_total_genes_present, 3),  # Percent total genes present
                num_unique_genes,  # Number of genes unique to group
                round(percent_unique_genes, 3),  # Percent unique genes
                num_genes_absent,  # Number of total genes absent from group
                round(percent_genes_absent, 3)  # Percent genes absent
            ])

    # Create DataFrame for output
    output_df = pd.DataFrame(results, columns=[
        'Group', 'Class', 'Total genomes in Group', 'Total genes in Group',
        'Percent total genes present', 'Number of genes unique to group', 'Percent unique genes',
        'Number of total genes absent from group', 'Percent genes absent'
    ])

    # Save output to file
    output_df.to_csv(f"{prefix}.unique_group_differences.tab", sep='\t', index=False)


def main(metadata_file, metadata2, presence_absence_file, prefix):
    # Load the data
    metadata, presence_absence = load_data(metadata_file, metadata2, presence_absence_file)

    # First output file - count/percentage of unique genes
    calculate_group_gene_metrics(metadata, presence_absence, prefix)

    # Instantiate the 2nd and 3rd output files - more detailed pair comparisons
    with open(f"{prefix}.pairwise_group_differences.tab", 'w') as out:
        with open(f"{prefix}.group_union_genes.tab", 'w') as file_unions:
            with open(f"{prefix}.group_unique_genes.tab", 'w') as file_uniques:
                
                # Write out headers for both files
                out.write("Class\tGroupA\tGroupB\tUnion\tIntersection\tAminusB\tBminusA\tSymDifference\tCardinalityA\tCardinalityB\n")
                file_unions.write("Gene\tClass\tLabel\n")
                file_uniques.write("Gene\tClass\tLabel\n")

                # Get the unique labels in the metadata column
                for metadata_column in list(metadata.columns.values):
                    labels = metadata[metadata_column].unique()

                    # Compute gene sets for each label
                    label_to_genes = {}
                    for label in labels:
                        genomes = metadata[metadata[metadata_column] == label].index.tolist()
                        genes = compute_gene_sets(presence_absence, genomes)
                        label_to_genes[label] = genes

                    # Compute pairwise union, difference, and symmetric difference
                    for i, label_a in enumerate(labels):
                        for j, label_b in enumerate(labels):
                            if i < j:  # To avoid redundancy (computing pairs twice)
                                genes_a = label_to_genes[label_a]
                                genes_b = label_to_genes[label_b]

                                union_genes = compute_union([genes_a, genes_b])
                                intersect_genes = compute_intersection([genes_a, genes_b])
                                adiffb_genes = compute_difference(genes_a, genes_b)
                                bdiffa_genes = compute_difference(genes_b, genes_a)
                                sym_diff_genes = compute_symmetric_difference(genes_a, genes_b)

                                # Write the main output line
                                out.write(f"{metadata_column}\t{label_a}\t{label_b}\t{len(union_genes)}\t{len(intersect_genes)}\t{len(adiffb_genes)}\t{len(bdiffa_genes)}\t{len(sym_diff_genes)}\t{len(genes_a)}\t{len(genes_b)}\n")

                    # Compute comparison of each label with the union of all other labels
                    for label in labels:
                        genomes_a = metadata[metadata[metadata_column] == label].index.tolist()
                        genes_a = label_to_genes[label]

                        # Create the union of all other labels except the current one
                        other_labels = [l for l in labels if l != label]
                        if len(other_labels) == 0: 
                            continue
                        else:
                            genomes_other = metadata[metadata[metadata_column].isin(other_labels)].index.tolist()
                            genes_other = compute_union([label_to_genes[l] for l in other_labels])

                            union_genes = compute_union([genes_a, genes_other])
                            intersect_genes = compute_intersection([genes_a, genes_other])
                            adiffb_genes = compute_difference(genes_a, genes_other)
                            bdiffa_genes = compute_difference(genes_other, genes_a)
                            sym_diff_genes = compute_symmetric_difference(genes_a, genes_other)

                            # Write the comparison of the current label with the union of all others
                            # This is slightly redundant in cases where there are only two categories,
                            #    but I'm opting to allow it, with an eye towards easier parsing/slicing downstream.
                            #    Change the >= to just > if you want to disallow this.
                            if len(other_labels) >= 1:
                                out.write(f"{metadata_column}\t{label}\tAll others in class\t{len(union_genes)}\t{len(intersect_genes)}\t{len(adiffb_genes)}\t{len(bdiffa_genes)}\t{len(sym_diff_genes)}\t{len(genes_a)}\t{len(genes_other)}\n")

                            # Write out the union genes into it's list file
                            for union_gene in union_genes:
                                file_unions.write(f"{union_gene}\t{metadata_column}\t{label}\n")
                            
                            # Write out the unique genes into our list file
                            for unique_gene in adiffb_genes:
                                file_uniques.write(f"{unique_gene}\t{metadata_column}\t{label}\n")


if __name__ == "__main__":
    # Parse incoming args
    args = parse_arguments()
    metadata_file = args.meta
    metadata_file2 = args.meta2
    presence_absence_file = args.matrix

    # Activate main loop
    main(metadata_file, metadata_file2, presence_absence_file, args.output)

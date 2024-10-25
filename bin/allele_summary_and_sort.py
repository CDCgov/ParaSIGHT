#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import re

def parse_subtype_family(subtype):
    if pd.isna(subtype) or subtype.strip() == '':
        return np.nan
    if re.match(r"^Assemblage", subtype):
        return ""
    match = re.match(r"^[A-Z]+[a-z]?", subtype)
    return match.group(0) if match else np.nan

def categorize_frequency(frequency):
    """
    Categorize allele frequency based on predefined thresholds.
    """
    if frequency >= 95:
        return "Extremely Common (100-95%)"
    elif frequency >= 85:
        return "Very Common (85-95%)"
    elif frequency >= 70:
        return "Somewhat Common (70-85%)"
    elif frequency >= 30:
        return "Common (30-70%)"
    elif frequency < 30 and frequency >= 15:
        return "Somewhat Uncommon (15-30%)"
    elif frequency < 15 and frequency >= 5:
        return "Very Uncommon (5-15%)"
    else:
        return "Extremely Uncommon (1-5%)"

# Set some constants / lookup tables
month_labels = {
        1: 'January', 2: 'February', 3: 'March', 4: 'April', 
        5: 'May', 6: 'June', 7: 'July', 8: 'August',
        9: 'September', 10: 'October', 11: 'November', 12: 'December',
        0: 'Unknown'
}
quarter_labels = {
        1: 1, 2: 1, 3: 1, 
        4: 2, 5: 2, 6: 2, 
        7: 3, 8: 3, 9: 3, 
        10: 4, 11: 4, 12: 4,
        0: "Unknown"
}
hhs_regions = {
    'AL': 4, 'AK': 10, 'AZ': 9, 'AR': 6, 'CA': 9, 'CO': 8, 'CT': 1, 'DE': 3,
    'FL': 4, 'GA': 4, 'HI': 9, 'ID': 10, 'IL': 5, 'IN': 5, 'IA': 7, 'KS': 7,
    'KY': 4, 'LA': 6, 'ME': 1, 'MD': 3, 'MA': 1, 'MI': 5, 'MN': 5, 'MS': 4,
    'MO': 7, 'MT': 8, 'NE': 7, 'NV': 9, 'NH': 1, 'NJ': 2, 'NM': 6, 'NY': 2,
    'NC': 4, 'ND': 8, 'OH': 5, 'OK': 6, 'OR': 10, 'PA': 3, 'RI': 1, 'SC': 4,
    'SD': 8, 'TN': 4, 'TX': 6, 'UT': 8, 'VT': 1, 'VA': 3, 'WA': 10, 'WV': 3,
    'WI': 5, 'WY': 8, 'DC': 3, 'NA': 'International'
}

# Main loop begins
parser = argparse.ArgumentParser(description="Helper script to merge metadata and outputs from alignment2alleles to CSV for CryptoNet's Alleles PowerBI dashboard.")
parser.add_argument("-m", "--meta", dest="meta", help="Input file path to metadata file")
parser.add_argument("-n", "--meta2", dest="meta2", help="(Optional) Second input file path to metadata file")
parser.add_argument("-a", "--alleles", dest="alleles", help="Input file path to `alleles_per_genome.tab` from alignment2alleles script")
parser.add_argument("-b", "--annots", dest="annots", help="(Optional) Input file path to `union_genes_annotation.tab` from crosstab_unique_genes script")
parser.add_argument("-g", "--genes", dest="genes", help="Input file path to `distribution.tab` from blastn2coregenome script")
parser.add_argument("-o", "--output", dest="file_out", help="Output file prefix (no extensions!)")
args = parser.parse_args()

# Load the allele and metadata files
allele_df = pd.read_csv(args.alleles, sep="\t")
metadata_df = pd.read_csv(args.meta)

# Merge additional metadata if provided
if args.meta2:
    metadata2_df = pd.read_csv(args.meta2)
    metadata_df = pd.concat([metadata_df, metadata2_df])
    metadata_df.drop(columns=['fasta'], inplace=True)

# Merge allele and metadata dataframes
allele_df = allele_df.merge(metadata_df, left_on='Genome', right_on='sample', how='inner')

# Load and merge genes information
genes_df = pd.read_csv(args.genes)
genes_df.drop(columns=['Total Genomes', 'Gene Length', 'No. In-group Genomes', 'In-group Mean', 'In-group Stdev', 'No. Across-group Genomes', 'Across-group Mean', 'Across-group Stdev', 'Failed QC'], inplace=True)
allele_df = allele_df.merge(genes_df, left_on='Gene', right_on='Gene', how='inner')
allele_df.drop(columns=['sample', 'fastq_1', 'fastq_2'], inplace=True)

# If annotations are provided, merge them too
if args.annots:
    annots_df = pd.read_csv(args.annots, sep="\t")
    annots_df['COG Category'] = annots_df['COG Category'].astype(str).fillna("Unannotated")
    annots_df['KEGG KO'] = annots_df['KEGG KO'].astype(str).fillna("Unannotated")
    annots_df.replace("NA", "", inplace=True)
    annots_df.replace("nan", "", inplace=True)
    annots_df.drop(columns=['Class', 'Label'], inplace=True)
    annots_df.drop_duplicates(inplace=True)

# Ensure the data types are correct
allele_df['year'] = allele_df['year'].astype(int).fillna('Unknown')
allele_df['month'] = allele_df['month'].astype(int).fillna('Unknown')

# Add in some additional data columns, extrapolated from existing input data
allele_df['monthName'] = allele_df["month"].map(month_labels)
allele_df['quarter'] = allele_df["month"].map(quarter_labels)
allele_df['hhsregion'] = allele_df['state'].map(hhs_regions, na_action='ignore').fillna("International")
allele_df['state'] = allele_df['state'].fillna(allele_df['country'])
allele_df['subtype_family'] = allele_df['subtype'].map(parse_subtype_family)
allele_df['Outbreak associated isolates'] = np.where(allele_df['casetype'].isin(["Outbreak", 'outbreak']), 1, 0)

# Sort the data by Year, Quarter, then Month, then State
df_sorted = allele_df.sort_values(by=['year', 'quarter', 'month', 'state'])

# Summarize the data
summary = df_sorted.groupby(['year', 'quarter', 'month', 'monthName', 'hhsregion', 'state', 'pathogen', 'species', 'subtype_family', 'subtype', 'Category', 'Gene', 'Allele' ]).agg({
    'Genome': 'nunique',
    'Outbreak associated isolates': 'sum'
}).reset_index()

# Calculate allele frequencies
allele_frequency_df = summary.groupby(['Gene', 'Allele'])['Genome'].nunique().reset_index()
allele_frequency_df.columns = ['Gene', 'Allele', 'Number of isolates']
total_genomes = summary['Genome'].nunique()

# Compute allele frequency as a percentage and categorize
allele_frequency_df['Frequency'] = np.round((allele_frequency_df['Number of isolates'] / total_genomes) * 100, 4)
allele_frequency_df['Allele Category'] = allele_frequency_df['Frequency'].apply(categorize_frequency)
allele_frequency_df.drop(columns=['Number of isolates'], inplace=True)
summary = summary.merge(allele_frequency_df, on=['Gene', 'Allele'], how='inner')

# Write first output files
if args.annots:
    summary = summary.merge(annots_df, left_on='Gene', right_on='Gene', how='inner')
    summary.columns = ['Year', 'Quarter', 'Month', 'MonthName', 'HHSregion', 'State', 'Pathogen', 'Species', 'Subtype family', 'Subtype', 'Category', 'Gene', 'Allele', 'Number of isolates', 'Outbreak associated isolates', 'Allele Frequency', 'Allele Category', 'COG Category', 'COG Class', "Gene Name", "Gene Description", 'KEGG KO', 'KEGG Module', 'KEGG Pathway', 'KEGG Path']
else:
     summary.columns = ['Year', 'Quarter', 'Month', 'MonthName', 'HHSregion', 'State', 'Pathogen', 'Species', 'Subtype family', 'Subtype', 'Category', 'Gene', 'Allele', 'Number of isolates', 'Outbreak associated isolates', 'Allele Frequency', 'Allele Category' ]
summary.to_csv(f"{args.file_out}.alleles_dashboard.csv", index=False)

# Print the output summary per genome file
summary2 = df_sorted.sort_values(by=['year', 'quarter', 'month', 'monthName', 'hhsregion', 'state', 'pathogen', 'species', 'subtype_family', 'subtype', 'Category', 'Gene', 'Allele', 'Genome' ])
summary2 = summary2.merge(allele_frequency_df, on=['Gene', 'Allele'], how='inner')
summary2 = summary2.merge(annots_df, left_on='Gene', right_on='Gene', how='inner')
summary2.to_csv(f"{args.file_out}.alleles_per_genome_dashboard.csv", index=False)
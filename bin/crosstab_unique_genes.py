#!/usr/bin/env python3

# Acknowledgement that some code in this script is modified from Github:
# Credit: https://github.com/rotheconrad/EggNog_Annotation_Plots/blob/main/scripts/annotation_bar_plot.py

import argparse
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from itertools import product


def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze gene content differences between groups of genomes.")
    parser.add_argument("-i", "--genes", required=True, help="Input table of unique genes plus their metadata file (in TAB format).")
    parser.add_argument("-a", "--annots", required=True, help="Input eggNOG emapper.annotations (in TAB format).")
    parser.add_argument("-k", "--keggs", required=True, help="Input KEGG annotations file (in TAB format).")
    parser.add_argument("-o", "--output", required=True, help="Output file prefix")
    return parser.parse_args()


# Function to read and merge data
def read_and_merge(gene_file, annotation_file):
    # Read the gene metadata file
    gene_df = pd.read_csv(gene_file, sep='\t')
    
    
    # Read the functional annotation file
    annotation_df = pd.read_csv(annotation_file, sep='\t', skiprows=5, header=None)
    annotation_df.columns = [ 'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                       'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 
                       'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                       'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                       'BiGG_Reaction', 'PFAMs' ]

    # Some minor data cleaning for the COGs annotations
    annotation_df['COG_category'].replace('-', "Not Assigned", inplace=True)   # Converts missing/unknown to S category ("Unknown Function")
    annotation_df['COG_category'].replace(r"[A-Z]{2,}", "Multiple COG", regex=True, inplace=True)   # Converts multi-COG to 'MultiCOG' category 

    # Merge the two tables on the 'Gene' and 'query_name' column
    merged_df = pd.merge(gene_df, annotation_df, left_on='Gene', right_on='query', how='inner')
    
    return merged_df

# Read in and parse the KEGG file into a lookup table
def parse_keggfile(keggfile):
    # {ko: [pathway, module, path]}
    keggs = {}
    with open(keggfile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            pathway = X[0]
            module = X[1]
            path = X[2]
            ko = X[3]
            keggs[ko] = [pathway, module, path]

    return keggs


# Read in and parse the `unique_genes` metadata file into lookup dictionary by gene name
def parse_metadata(metadata):
    md = {}
    with open(metadata, 'r') as file:
        for line in file:
            # Skip header
            if line.startswith("Gene\t"): 
                continue
            
            # Split the line by tab
            X = line.rstrip().split('\t')
            
            # Assign values for sequence name, class name, and label name
            sn = X[0]   # sequence name (Gene)
            gn = X[1]   # class name (Class)
            ln = X[2]   # label name (Label)
 
            # Assign class name as key and label name as value in the dictionary
            md[sn] = [ gn, ln ]

    return md


def parse_annotations(annotations, md, keggs):
    # Instantiate ano dictionary
    ano = {
            'Sequence Name': [],
            'Class Name': [],
            'Label Name': [],
            'COG Category': [],
            'COG Class': [],
            'Gene Name': [],
            'Gene Description': [],
            'KEGG KO': [],
            'KEGG Pathway': [],
            'KEGG Module': [],
            'KEGG Path': [],
            }
    
    # Assign COG category
    COGclass = {
            'X': 'Mobile', 'C': 'Metabolism 1', 'G': 'Metabolism 1',
            'E': 'Metabolism 1', 'F': 'Metabolism 1', 'H': 'Metabolism 1',
            'I': 'Metabolism 1', 'P': 'Metabolism 2', 'Q': 'Metabolism 2',
            'J': 'Ribosomal', 'A': 'Information', 'K': 'Information',
            'L': 'Information', 'B': 'Information', 'D': 'Cellular',
            'Y': 'Cellular', 'V': 'Cellular', 'T': 'Cellular',
            'M': 'Cellular', 'N': 'Cellular', 'Z': 'Cellular',
            'W': 'Cellular', 'U': 'Cellular', 'O': 'Cellular',
            'S': 'Conserved Hypothetical', 'R': 'Conserved Hypothetical',
            '-': 'Hypothetical'
            }

    # keywords to assign category for EggNog
    mobile = [
                'transposase', 'phage', 'integrase', 'viral', 'plasmid',
                'integron', 'transposon'
                ]

    with open(annotations, 'r') as file:
        for line in file:
            if line.startswith('#'): 
                continue
            X = line.rstrip().split('\t')
            name = X[0]                         # representitive predicted gene name
            if name not in md:
                continue
            else:
                group_name = md[name][0]
                label_name = md[name][1]
                gene = X[8]                     # annotation short gene name
                if gene == '-': gene = 'NA'
                desc = X[7]                     # annotation long gene name (description)
                if desc == '-': desc = 'NA'
                cog = X[6][0]                   # select only first letter
                ko = X[11].split(',')[0]
                if ko == '-': ko = 'Hypothetical'
                else: ko = ko.split(':')[1]
                if keggs and ko != 'Hypothetical':
                    d = keggs.get(ko, ['NA', 'NA', 'NA'])
                    pathway, module, path = d[0], d[1], d[2]
                else:
                    pathway, module, path = 'NA', 'NA', 'NA'

                if any(mbl in desc.lower() for mbl in mobile):
                    cog, cat = 'X', 'Mobile'
                elif cog == '-': cog, cat = 'NA', 'Hypothetical'
                else: cat = COGclass.get(cog, 'Other')

            # Populate the gene's data in the ano dictionary
            ano['Sequence Name'].append(name)
            ano['Class Name'].append(group_name)
            ano['Label Name'].append(label_name)
            ano['COG Category'].append(cat)
            ano['COG Class'].append(cog)
            ano['Gene Name'].append(gene)
            ano['Gene Description'].append(desc)
            ano['KEGG KO'].append(ko)
            ano['KEGG Pathway'].append(pathway)
            ano['KEGG Module'].append(module)
            ano['KEGG Path'].append(path)

    df = pd.DataFrame(ano)
    return df


# Function to create crosstabulation and plots
def crosstabulate_and_plot(merged_df, prefix):
    classes = merged_df['Class'].unique()
    
    for class_name in classes:
        class_df = merged_df[merged_df['Class'] == class_name]
        
        # Crosstabulate 'COG_category' with 'Label'
        crosstab_df = pd.crosstab(class_df['COG_category'], class_df['Label'])
        
        # Plot pie chart of COG category
        #plot_pie_chart(class_df, class_name, prefix)


# Function to plot pie chart for COG categories
# def plot_pie_chart(class_df, class_name, prefix):
#     # Count the number of occurrences for each 'COG_category'
#     cog_counts = class_df['COG_category'].value_counts()

#     # Keep the top 10 COG categories, and group the rest as "Other"
#     top_10_cogs = cog_counts.nlargest(10)
#     other_cogs = cog_counts.iloc[10:].sum()  # Sum all remaining COG categories
    
#     # Add the "Other" category if there are any remaining COGs
#     if other_cogs > 0:
#         top_10_cogs["All Other COGs"] = other_cogs
    
#     # Create pie chart
#     plt.figure(figsize=(8, 8))
#     plt.pie(top_10_cogs, labels=top_10_cogs.index, autopct='%1.1f%%', startangle=90, colors=plt.cm.Paired.colors)
#     plt.title(f"COG Categories Distribution for Class: {class_name}")
#     plt.tight_layout()
#     plt.savefig(f"{prefix}_{class_name}_cog-pie_mqc.png")
#     plt.close()


# Main function to run the entire pipeline
def main(gene_file, annotation_file, kegg_file, prefix):
    
    # Read in the various data files
    keggs = parse_keggfile(kegg_file) if kegg_file else None
    md = parse_metadata(gene_file)
    annots = parse_annotations(annotation_file, md, keggs)

    # Merge the annotations with the metadata
    gene_df = pd.read_csv(gene_file, sep='\t')
    blend_df = pd.merge(gene_df, annots, left_on=['Gene'], right_on=['Sequence Name'], how='left')
    #blend_df.replace(np.nan, "NA", inplace=True)
    blend_df = blend_df.drop(columns=['Sequence Name', 'Class Name', 'Label Name'])
    blend_df.to_csv(f"{prefix}_gene_annotations.tab", sep='\t', index=False)

    # Crosstabulation and plotting for unique genes
    merged_df = read_and_merge(gene_file, annotation_file)
    crosstabulate_and_plot(merged_df, prefix)


# Example usage
if __name__ == '__main__':
    # Parse incoming args
    args = parse_arguments()
    gene_file = args.genes 
    annotation_file = args.annots
    kegg_file = args.keggs
    prefix = args.output
    
    main(gene_file, annotation_file, kegg_file, prefix)

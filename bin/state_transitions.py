#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from itertools import combinations

# Function to calculate Jaccard similarity with pseudocounts
def jaccard_similarity_with_pseudocount(set1, set2):
    # Add pseudocount of 1 to both sets
    intersection = len(set1.intersection(set2)) + 1  # Adding pseudocount to intersection
    union = len(set1.union(set2)) + 1  # Adding pseudocount to union
    return intersection / union

# Function to compute pairwise Jaccard similarity matrix for each category
def calculate_jaccard(df, category, allele_column='Allele'):
    # Get unique values for the category
    unique_values = df[category].unique()

    # Error handling for cases with only one unique value
    if len(unique_values) < 2:
        print(f"Skipping category '{category}' because it only has one unique value.")
        return None
    
    # Initialize an empty matrix (pairwise comparison matrix)
    matrix = pd.DataFrame(index=unique_values, columns=unique_values)
    
    # Calculate Jaccard similarity for each pair with pseudocounts
    for val1, val2 in combinations(unique_values, 2):
        set1 = set(df[df[category] == val1][allele_column])
        set2 = set(df[df[category] == val2][allele_column])
        
        similarity = jaccard_similarity_with_pseudocount(set1, set2)
        
        matrix.at[val1, val2] = similarity
        matrix.at[val2, val1] = similarity  # Symmetric matrix

    # Fill diagonal with 1 (self-comparison)
    for val in unique_values:
        matrix.at[val, val] = 1

    # Normalize each row to sum to 1
    matrix = np.round(matrix.div(matrix.sum(axis=1), axis=0), 4)

    return matrix


def main(input_csv, output_csv):
    # Read in the CSV file
    df = pd.read_csv(input_csv)

    # List of categories for which you want to compute Jaccard similarity
    categories = ['pathogen', 'species', 'subtype', 'country', 'state', 'year', 
                'host', 'source', 'casetype', 'monthName', 'quarter', 
                'hhsregion', 'subtype_family']

    # Calculate and output Jaccard similarity for each category
    output_matrices = {}

    for category in categories:
        print(f"Calculating Jaccard similarity for {category}...")
        jaccard_matrix = calculate_jaccard(df, category)

        if jaccard_matrix is not None:  # Skip saving if there is no matrix
            output_matrices[category] = jaccard_matrix
            
            # Save the matrix to a CSV file
            output_file = f"{category}.state_transitions.csv"
            jaccard_matrix.to_csv(output_file, index=True)
            print(f"Saved Jaccard similarity matrix for {category} to {output_file}")

if __name__ == "__main__":
    # Argument parser for input and output files
    parser = argparse.ArgumentParser(description="Calculate a state transition matrix for discrete categories and output as a table.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    # Run the main function
    main(args.input, args.output)

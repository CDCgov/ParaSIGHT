#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np

def calculate_allele_frequencies(df, grouping_columns):
    """
    Calculate allele frequencies by grouping by Allele and the specified categories.

    Parameters:
        df (pd.DataFrame): The input DataFrame.
        grouping_columns (list): List of columns to group by and calculate frequencies.

    Returns:
        pd.DataFrame: A DataFrame where each row is an allele and each column is the category with allele frequencies.
    """
    # Initialize a dictionary to hold frequency data
    allele_freq_dict = {}

    # Loop through each grouping category
    for category in grouping_columns:
        # Group by Allele and the given category, and calculate frequency within each group
        group_df = df.groupby(['Allele', category]).size().reset_index(name=f'{category} Frequency')
        # Normalize frequency for each Allele
        total_counts = df.groupby(['Allele']).size().reset_index(name='Total Count')
        group_df = pd.merge(group_df, total_counts, on=['Allele'])
        group_df[f'{category} Frequency'] = np.round(group_df[f'{category} Frequency'] / group_df['Total Count'], 4)
        
        # Pivot the data to get categories as columns
        pivot_df = group_df.pivot_table(index=['Allele'], columns=category, values=f'{category} Frequency', fill_value=0)
        
        # Flatten the MultiIndex columns
        pivot_df.columns = [f'{category}_{col}' for col in pivot_df.columns]
        
        # Add to the dictionary
        allele_freq_dict[category] = pivot_df

    # Concatenate all frequency tables together on Gene and Allele
    final_df = pd.concat(allele_freq_dict.values(), axis=1)

    return final_df.reset_index()

def main(input_csv, output_csv):
    # Read in the CSV file
    df = pd.read_csv(input_csv)
    
    # Define the columns for which we need to calculate allele frequencies
    grouping_columns = ['pathogen', 'species', 'subtype', 'country', 'state', 'year', 
                        'host', 'source', 'casetype', 'monthName', 'quarter', 'hhsregion', 
                        'subtype_family']

    # Calculate allele frequencies
    allele_freq_df = calculate_allele_frequencies(df, grouping_columns)

    # Write the result to a CSV file
    allele_freq_df.to_csv(output_csv, index=False)

    print(f"Allele frequency table written to {output_csv}")

if __name__ == "__main__":
    # Argument parser for input and output files
    parser = argparse.ArgumentParser(description="Calculate allele frequencies for each category and output as a table.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    # Run the main function
    main(args.input, args.output)

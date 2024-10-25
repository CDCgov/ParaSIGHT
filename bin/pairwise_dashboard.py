#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from datetime import datetime, date

hhs_regions = {
    'AL': 4, 'AK': 10, 'AZ': 9, 'AR': 6, 'CA': 9, 'CO': 8, 'CT': 1, 'DE': 3,
    'FL': 4, 'GA': 4, 'HI': 9, 'ID': 10, 'IL': 5, 'IN': 5, 'IA': 7, 'KS': 7,
    'KY': 4, 'LA': 6, 'ME': 1, 'MD': 3, 'MA': 1, 'MI': 5, 'MN': 5, 'MS': 4,
    'MO': 7, 'MT': 8, 'NE': 7, 'NV': 9, 'NH': 1, 'NJ': 2, 'NM': 6, 'NY': 2,
    'NC': 4, 'ND': 8, 'OH': 5, 'OK': 6, 'OR': 10, 'PA': 3, 'RI': 1, 'SC': 4,
    'SD': 8, 'TN': 4, 'TX': 6, 'UT': 8, 'VT': 1, 'VA': 3, 'WA': 10, 'WV': 3,
    'WI': 5, 'WY': 8, 'DC': 3, 'NA': 'International'
}

def compute_time_comparison(date1, date2):
    """
    Computes the time comparison between two genome samples based on their dates using a rolling window.
    """
    delta = abs((date2 - date1).days)
    
    if delta <= 45:
        return "Approx. 30-45 days (same month)"
    elif delta <= 120:
        return "Less than 120 days"
    elif delta <= 365:
        return "Within 12 months"
    else:
        return "Longer than 12 months apart"

def compute_geo_comparison(state1, state2, hhs_region1, hhs_region2):
    """
    Computes the geographic comparison between two genome samples.
    """
    if state1 == state2:
        return "Same state"
    elif hhs_region1 == hhs_region2:
        return "Multistate - same HHS Region"
    else:
        return "Multistate - Widespread"

def compute_ani_category(ani):
    """
    Categorizes the ANI value into predefined bins.
    """
    if ani >= 99:
        return ">99"
    elif 99 > ani >= 95:
        return "99-95"
    elif 95 > ani >= 90:
        return "95-90"
    elif 90 > ani >= 85:
        return "90-85"
    elif 85 > ani >= 80:
        return "85-80"
    else:
        return "<80"

def compute_case_comparison(ani, mlst_jaccard, allele_jaccard):
    """
    Determines the case comparison based on ANI, MLST Jaccard, and Allele Jaccard values.
    """
    # The values below are empirically chosen and set to account for small rounding error
    if ani >= 99.850 and mlst_jaccard >= 0.9850 and allele_jaccard >= 0.9450:
        return "Strongly linked"
    elif ani >= 99.6850 and mlst_jaccard >= 0.9450 and allele_jaccard >= 0.8950:
        return "Possibly linked (relaxed linkage)"
    else:
        return "Unlikely to be linked"

def main():
    parser = argparse.ArgumentParser(description="Join genome comparisons with metadata and calculate labels.")
    parser.add_argument("-c", "--comparisons", required=True, help="Input file path to genome pairwise comparison table (csv)")
    parser.add_argument("-m", "--meta", required=True, help="Input file path to metadata table (csv)")
    parser.add_argument("-n", "--meta2", dest="meta2", help="(Optional) Second input file path to metadata file (csv)")
    parser.add_argument("-o", "--output", required=True, help="Output file path for the labeled CSV file")
    args = parser.parse_args()

    # Load the pairwise comparison table and metadata
    comp_df = pd.read_csv(args.comparisons)
    meta_df = pd.read_csv(args.meta)
    meta_df.drop(columns=['fastq_1', 'fastq_2'], inplace=True)

    # Merge additional metadata if provided
    if args.meta2:
        meta2_df = pd.read_csv(args.meta2)
        meta_df = pd.concat([meta_df, meta2_df])
        meta_df.drop(columns=['fasta'], inplace=True)

    # Parse metadata date columns into datetime objects --> use 01/01/1900 for missing/unknown dates
    meta_df['year'] = meta_df['year'].fillna('1900').astype(str)
    meta_df['month'] = meta_df['month'].fillna('0').astype(str)
    meta_df['date'] = pd.to_datetime(meta_df['year'].astype(str) + '-' + meta_df['month'].astype(str), format="%Y-%m")
    meta_df['hhsregion'] = meta_df['state'].map(hhs_regions, na_action='ignore').fillna("International")
    meta_df['state'] = meta_df['state'].fillna(meta_df['country'])
    
    # Join metadata on both Genome1 and Genome2
    merged_df = comp_df.merge(meta_df, left_on='Genome1', right_on='sample', suffixes=('_1', '_2'))
    merged_df = merged_df.merge(meta_df, left_on='Genome2', right_on='sample', suffixes=('_1', '_2'))
    
    print("\n\n")
    print(merged_df)
    merged_df.to_csv("merged.csv", sep="\t", index=False)
    print("\n\n")

    # Compute the comparison labels
    merged_df['Time Comparison'] = merged_df.apply(lambda row: compute_time_comparison(
        row['date_1'], row['date_2']), axis=1)

    merged_df['Geo Comparison'] = merged_df.apply(lambda row: compute_geo_comparison(
        row['state_1'], row['state_2'], row['hhsregion_1'], row['hhsregion_2']), axis=1)

    merged_df['ANI Category'] = merged_df['ANI'].apply(compute_ani_category)

    merged_df['Case Comparison'] = merged_df.apply(lambda row: compute_case_comparison(
        row['ANI'], row['Gene Jaccard'], row['Allele Jaccard']), axis=1)

    # Select the relevant columns for output
    output_columns = [
        'Genome1', 'Genome2', 'pathogen_1', 'species_1', 'subtype_1', 'hhsregion_1', 'state_1', 'year_1', 'month_1',
        'host_1', 'source_1', 'casetype_1', 'pathogen_2', 'species_2', 'subtype_2', 'hhsregion_2', 'state_2', 'year_2',
        'month_2', 'host_2', 'source_2', 'casetype_2', 'ANI', 'Gene Difference', 'Gene Jaccard', 'Allele Jaccard',
        'Time Comparison', 'Geo Comparison', 'ANI Category', 'Case Comparison'
    ]
    output_df = merged_df[output_columns]

    # Write the output to a CSV file
    output_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()

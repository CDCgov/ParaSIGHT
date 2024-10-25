#!/usr/bin/env python3

import csv
import sqlite3
import argparse

# Helper function to handle empty strings by converting them to None
def handle_null(value):
    return value if value != '' else None

# Function to insert or get Location ID
def get_or_insert_location(cursor, country, state):
    country = handle_null(country)
    state = handle_null(state)
    
    # Check if the location already exists
    cursor.execute("SELECT location_id FROM Location WHERE country = ? AND state = ?", (country, state))
    result = cursor.fetchone()
    
    # If found, return the location_id
    if result:
        return result[0]
    
    # Otherwise, insert a new location
    cursor.execute("INSERT INTO Location (country, state) VALUES (?, ?)", (country, state))
    return cursor.lastrowid

# Function to insert or get Host ID
def get_or_insert_host(cursor, host):
    host = handle_null(host)
    
    # Check if the host already exists
    cursor.execute("SELECT host_id FROM Host WHERE host = ?", (host,))
    result = cursor.fetchone()
    
    # If found, return the host_id
    if result:
        return result[0]
    
    # Otherwise, insert a new host
    cursor.execute("INSERT INTO Host (host) VALUES (?)", (host,))
    return cursor.lastrowid

# Function to insert data into the Genome_Metadata table
def insert_genome_metadata(cursor, data):
    # Prepare the data, handling empty strings as NULL
    csid = handle_null(data['csid'])
    cuid = handle_null(data['cuid'])
    original_id = handle_null(data['original_id'])
    mepi_id = handle_null(data['mepi_id'])
    pathogen = handle_null(data['pathogen'])
    species = handle_null(data['species'])
    subtype = handle_null(data['subtype'])
    source = handle_null(data['source'])
    casetype = handle_null(data['casetype'])
    
    # Get or insert Location ID
    location_id = get_or_insert_location(cursor, data['country'], data['state'])
    
    # Get or insert Host ID
    host_id = get_or_insert_host(cursor, data['host'])
    
    # Insert into Genome_Metadata table
    cursor.execute("""
        INSERT INTO Genome_Metadata (
            csid, cuid, original_id, mepi_id, pathogen, species, subtype, location_id, host_id, source, casetype
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (csid, cuid, original_id, mepi_id, pathogen, species, subtype, location_id, host_id, source, casetype))

# Main function to process the CSV and upload data
def process_csv_and_upload(csv_file, db_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Open the CSV file and process each row
    with open(csv_file, newline='') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            insert_genome_metadata(cursor, row)

    # Commit the changes and close the connection
    conn.commit()
    conn.close()
    print(f"Data from {csv_file} has been successfully uploaded to {db_file}.")

# Function to handle command-line arguments using argparse
def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Upload genome metadata from a CSV file into a SQLite database.")
    
    # Add arguments
    parser.add_argument("-i", "--input", dest="csv_file", required=True, help="The path to the CSV file containing genome metadata")
    parser.add_argument("-d", "--db", dest="db_file", required=True, help="The path to the SQLite database file")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function to process and upload CSV data
    process_csv_and_upload(args.csv_file, args.db_file)

# Entry point
if __name__ == '__main__':
    main()



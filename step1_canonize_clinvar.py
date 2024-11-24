import pandas as pd
import requests
from Bio.Seq import Seq
import sys
import os
import time

def fetch_canonical_cds_from_transcript(transcript_id):
    transcript_id = transcript_id.split('.')[0]  # Remove version number if present
    server = "https://rest.ensembl.org"
    endpoint = f"/sequence/id/{transcript_id}"
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(server + endpoint, headers=headers, params={"type": "cds"})
    if response.ok:
        sequence_data = response.json()
        return sequence_data["seq"]

    raise RuntimeError(f"Failed to fetch CDS for transcript {transcript_id}. Response: {response.text}")

def fetch_transcript_id(refseq_id, gene_symbol):
    try:
        server = "https://rest.ensembl.org"
        endpoint = f"/xrefs/id/{refseq_id}"
        headers = {"Content-Type": "application/json"}
        
        response = requests.get(server + endpoint, headers=headers)
        if response.ok:
            xrefs_data = response.json()
            for entry in xrefs_data:
                if entry.get("dbname") == "Ensembl_Transcript":
                    return entry.get("id")
    except Exception as e:
        print(f"Error in Ensembl lookup: {e}")

    try:
        server = "https://rest.ensembl.org"
        endpoint = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=true"
        headers = {"Content-Type": "application/json"}
        
        response = requests.get(server + endpoint, headers=headers)
        if response.ok:
            gene_data = response.json()
            return gene_data.get("canonical_transcript")
    except Exception as e:
        print(f"Error in gene lookup: {e}")

    raise RuntimeError(f"Could not resolve transcript ID for {refseq_id} or {gene_symbol}.")

def apply_mutation(cds_sequence, mutation_info):
    mutation_details = mutation_info.lstrip('c.')
    position, change = mutation_details[:-3], mutation_details[-3:]
    position = int(position) - 1  # Convert to 0-based index
    ref, alt = change.split('>')

    if cds_sequence[position] != ref:
        raise ValueError(f"Reference base mismatch at position {position + 1}: "
                         f"expected {ref}, found {cds_sequence[position]}.")
    mutated_cds = list(cds_sequence)
    mutated_cds[position] = alt
    return ''.join(mutated_cds)

def translate_cds(cds_sequence):
    cds = Seq(cds_sequence)
    protein = cds.translate(to_stop=True)
    return str(protein)

def find_protein_mutation(original_protein, mutated_protein):
    for i, (original, mutated) in enumerate(zip(original_protein, mutated_protein), start=1):
        if original != mutated:
            return f"{original}{i}{mutated}"
    return "No mutation found"

def process_mutation(mutation_info):
    """
    Processes a single mutation and returns the canonical mutation.
    """
    try:
        refseq_id = mutation_info.split('(')[0]
        gene_symbol = mutation_info.split('(')[1].split(')')[0]
        mutation = mutation_info.split(':')[1].split(' ')[0]

        transcript_id = fetch_transcript_id(refseq_id, gene_symbol)
        cds_sequence = fetch_canonical_cds_from_transcript(transcript_id)
        original_protein = translate_cds(cds_sequence)

        mutated_cds = apply_mutation(cds_sequence, mutation)
        mutated_protein = translate_cds(mutated_cds)

        return find_protein_mutation(original_protein, mutated_protein)
    except Exception as e:
        return f"Error: {e}"

def main():
    # Check if the input file is provided
    if len(sys.argv) < 2:
        print("Usage: python process_clinvar.py <clinvar_file>")
        sys.exit(1)

    # Get the input file from command-line arguments
    input_file = sys.argv[1]

    # Ensure the input file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    # Load the CSV file
    df = pd.read_csv(input_file)

    # Add a column for results and initialize with NaN
    df['Canonic mutation'] = pd.NA

    # Total rows for progress tracking
    total_rows = len(df)

    # Output file path
    output_file = os.path.join(os.getcwd(), "canonized_clinvar.csv")

    # Process each mutation and save the entire DataFrame after every query
    for idx in range(total_rows):
        print(f"Processing {idx + 1}/{total_rows}")
        mutation_info = df.loc[idx, 'Name']
        df.loc[idx, 'Canonic mutation'] = process_mutation(mutation_info)

        # Save the entire DataFrame to ensure progress is written
        df.to_csv(output_file, index=False)

    print(f"Final DataFrame saved to {output_file}")

if __name__ == "__main__":
    main()

import pandas as pd
from tqdm import tqdm

# File paths
clinvar_file = "canonized_clinvar.csv"
alphamissense_file = "/home/lucas/AlphaMissense_Benchmark/AlphaMissense_hg38.tsv"

# Load the dataframes
clinvar_df = pd.read_csv(clinvar_file)
alphamissense_df = pd.read_csv(alphamissense_file, sep='\t', skiprows=3)  # Skip the first 3 lines of AlphaMissense file

# Filter AlphaMissense to include only rows with chromosome 1 and limit rows for testing
#alphamissense_df = alphamissense_df[alphamissense_df['#CHROM'] == 'chr1']
#alphamissense_df = alphamissense_df[:1000000]

# Normalize columns to ensure consistent matching
clinvar_df['GRCh38Chromosome'] = clinvar_df['GRCh38Chromosome'].astype(str).str.strip()
alphamissense_df['#CHROM'] = alphamissense_df['#CHROM'].str.replace("chr", "").str.strip()

# Handle non-numeric GRCh38Location: keep only numeric rows, skipping invalid ones
clinvar_df['GRCh38Location'] = clinvar_df['GRCh38Location'].str.split('-').str[0].str.strip()  # Handle ranges
clinvar_df = clinvar_df[clinvar_df['GRCh38Location'].str.match(r'^\d+$', na=False)]  # Keep only numeric
clinvar_df['GRCh38Location'] = pd.to_numeric(clinvar_df['GRCh38Location'], errors='coerce')  # Convert to numbers

# Ensure protein variant columns are normalized
clinvar_df['Canonic mutation'] = clinvar_df['Canonic mutation'].str.strip().str.upper()
alphamissense_df['protein_variant'] = alphamissense_df['protein_variant'].str.strip().str.upper()

# Initialize new columns in the ClinVar dataframe
clinvar_df['am_pathogenicity'] = None
clinvar_df['am_class'] = None

# Iterate through the filtered AlphaMissense file using a progress bar
for _, am_row in tqdm(alphamissense_df.iterrows(), total=alphamissense_df.shape[0], desc="Processing AlphaMissense"):
    chrom = am_row['#CHROM']
    pos = am_row['POS']
    protein_variant = am_row['protein_variant']

    # Match in the ClinVar dataframe
    matches = clinvar_df[
        (clinvar_df['GRCh38Chromosome'] == chrom) &  # Match chromosome
        (clinvar_df['GRCh38Location'] == pos) &  # Match position
        (clinvar_df['Canonic mutation'] == protein_variant)  # Match protein variant
    ]

    # If matches are found, update the new columns
    if not matches.empty:
        indices = matches.index
        clinvar_df.loc[indices, 'am_pathogenicity'] = am_row['am_pathogenicity']
        clinvar_df.loc[indices, 'am_class'] = am_row['am_class']

# Save the updated ClinVar dataframe
output_file = "updated_clinvar_alphamiss.csv"
clinvar_df.to_csv(output_file, index=False)

print(f"Updated file saved to {output_file}")

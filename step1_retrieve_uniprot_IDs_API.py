import pandas as pd
import re
import numpy as np

from unipressed import UniprotkbClient, UniparcClient
import time
from unipressed import IdMappingClient
import pandas as pd
import requests


import csv
import subprocess
import os
import subprocess
from io import StringIO

from concurrent.futures import ThreadPoolExecutor

########
### Params
#######

run_grep = False


clinvar = pd.read_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_result.txt', sep='\t')
significance = clinvar['Clinical significance (Last reviewed)']
unique_genes = list(clinvar['Gene(s)'].unique())

final_genes = []
for i in range(len(unique_genes)):
    tmp_seq = unique_genes[i].split('|')
    for s in tmp_seq:
        final_genes.append(s)

unique_genes = list(set(final_genes))

def retrieve_uniprot_ids(gene_IDs, csv_path, cont=False):
    columns = ["from", "to"]

    final_data = pd.DataFrame(columns=columns)
    ranIDs = []
    uniprot_ids = []

    if cont == False:
        for gene in gene_IDs:


            print(gene)

            request = IdMappingClient.submit(
                source="GeneCards", dest="UniProtKB", ids={gene}
            )
            time.sleep(5)

            tmp=pd.DataFrame(list(request.each_result()))
            if len(tmp) != 0:

                ranIDs.append(gene)

                uniprot_ids.append(tmp['to'][0])

                # Create a DataFrame for the current iteration
                iteration_df = pd.DataFrame({
                    'from': ranIDs,
                    'to': uniprot_ids
                })

                # Save the DataFrame to a CSV file
                iteration_df.to_csv(csv_path, index=False)

            else:

                ranIDs.append(gene)

                uniprot_ids.append('Not found')

                # Create a DataFrame for the current iteration
                iteration_df = pd.DataFrame({
                    'from': ranIDs,
                    'to': uniprot_ids
                })

                # Save the DataFrame to a CSV file
                iteration_df.to_csv(csv_path, index=False)

    if cont == True:
        final_data = pd.read_csv(csv_path)
        #final_data.drop('Unnamed: 0', axis=1)
        ranIDs = list(final_data['from'])
        uniprot_ids = list(final_data['to'])

        gene_ID_new = [x for x in gene_IDs if x not in list(final_data['from'])]
        gene_IDs = gene_ID_new


        for gene in gene_IDs:


            print(gene)

            request = IdMappingClient.submit(
                source="GeneCards", dest="UniProtKB", ids={gene}
            )
            time.sleep(5)

            tmp=pd.DataFrame(list(request.each_result()))
            if len(tmp) != 0:

                ranIDs.append(gene)

                uniprot_ids.append(tmp['to'][0])

                # Create a DataFrame for the current iteration
                iteration_df = pd.DataFrame({
                    'from': ranIDs,
                    'to': uniprot_ids
                })

                # Save the DataFrame to a CSV file
                iteration_df.to_csv(csv_path, index=False)

            else:

                ranIDs.append(gene)

                uniprot_ids.append('Not found')

                # Create a DataFrame for the current iteration
                iteration_df = pd.DataFrame({
                    'from': ranIDs,
                    'to': uniprot_ids
                })

                # Save the DataFrame to a CSV file
                iteration_df.to_csv(csv_path, index=False)


    return 'Finished'

### Retrieve Uniprot IDs
cata_IDs = retrieve_uniprot_ids(unique_genes, csv_path="/home/lucas/cataract/cataract_AlphaMiss/uniprotIDs.csv", cont=True)

###################################################
## Create table with - Gene - Uniprot - Missense
###################################################

cata_IDs = pd.read_csv('/home/lucas/cataract/cataract_AlphaMiss/uniprotIDs.csv')

to_take = []
for i in range(len(clinvar)):
    if len(clinvar['Gene(s)'][i].split('|')) == 1:
        to_take.append(i)


clinvar_mono = clinvar[clinvar.index.isin(to_take)]

uniprot_column = []

for gene in clinvar_mono['Gene(s)']:
    for i in range(len(cata_IDs)):
        if cata_IDs['from'][i] == gene:
            uniprot_column.append(cata_IDs['to'][i])



def process_and_fgrep_batch(csv_file, tsv_file, output_dir, batch_size=80):
    unique_values = set()

    with open(csv_file, 'r') as csv_file:
        # Create a CSV reader
        csv_reader = csv.reader(csv_file)

        # Skip the header
        next(csv_reader, None)

        # Collect unique values from the second column
        for row in csv_reader:
            unique_values.add(row[1])

    # Split unique values into batches
    value_batches = [list(unique_values)[i:i + batch_size] for i in range(0, len(unique_values), batch_size)]

    # Process each batch in parallel
    with ThreadPoolExecutor() as executor:
        for batch in value_batches:
            # Use executor.map to parallelize processing for each value in the batch
            results = list(executor.map(run_fgrep, batch, [tsv_file] * len(batch), [output_dir] * len(batch)))

            # Process the results if needed

def run_fgrep(value, tsv_file, output_dir):
    try:
        fgrep_result = subprocess.check_output(['fgrep', value, tsv_file], universal_newlines=True)
    except subprocess.CalledProcessError as e:
        fgrep_result = ""

    if fgrep_result:
        output_file_path = os.path.join(output_dir, f"{value}_output.tsv")
        with open(output_file_path, 'w') as output_file:
            output_file.write(fgrep_result)
        print(f"Processing completed for: {value}")
    else:
        print(f"No matches found for: {value}")

# Example usage:
csv_file_path = '/home/lucas/cataract/cataract_AlphaMiss/uniprotIDs.csv'
tsv_file_path = '/home/lucas/alphamiss/AlphaMissense_aa_substitutions.tsv'
output_directory = '/home/lucas/cataract/cataract_AlphaMiss/database/'

if run_grep == True:
    process_and_fgrep_batch(csv_file_path, tsv_file_path, output_directory)

###########################
##### 
###########################

##############
### Remove fields with nan in mutations
################


clinvar_mono['uniprot'] = uniprot_column
clinvar_mono['Alpha_miss'] = None
clinvar_mono['Alpha_score'] = None
clinvar_mono['aa_references'] = None


clinvar_mono = clinvar_mono.dropna(subset=['Protein change'])

all_pheno = []
all_scores = []
all_references = []
uniprot_entries = list(clinvar_mono['uniprot'])
protein_changes = list(clinvar_mono['Protein change'])

verify = []
for i in range(len(uniprot_entries)):


    gene = uniprot_entries[i]
    mutations = protein_changes[i].split(', ')

    try:
        file_path = '/home/lucas/cataract/cataract_AlphaMiss/database/'+gene+'_output.tsv'
        alpha = pd.read_csv(file_path, sep='\t', header=None)

        tmp_pheno = []
        tmp_score = []
        tmp_reference_aa = []

        scores =  list(alpha[2])
        phenotypes = list(alpha[3])


        for clin_mut in range(len(mutations)):
            for alpha_mut in range(len(alpha)):
                if mutations[clin_mut][1:] == alpha[1][alpha_mut][1:]:

                    tmp_pheno.append(phenotypes[alpha_mut])
                    tmp_score.append(scores[alpha_mut])

                    if mutations[clin_mut] == alpha[1][alpha_mut]:
                        tmp_reference_aa.append('same reference')
                    else:
                        tmp_reference_aa.append('different reference: '+alpha[1][alpha_mut][0])


                if len(tmp_pheno) != len(tmp_score):
                   print('difference inside'), print(gene), print(mutations[clin_mut])


        all_pheno.append(tmp_pheno)
        all_scores.append(tmp_score)
        all_references.append(tmp_reference_aa)

        if len(all_pheno) != len(all_scores):
           break



    except FileNotFoundError as e:
        error_message = f"Error: File not found at {file_path}"
        all_pheno.append(error_message)
        all_scores.append(error_message)
        all_references.append(error_message)

    except pd.errors.EmptyDataError as e:
        error_message = f"Error: The file at {file_path} is empty or contains only whitespace."
        all_pheno.append(error_message)
        all_scores.append(error_message)
        all_references.append(error_message)


    except pd.errors.ParserError as e:
        error_message = f"Error: Unable to parse the file at {file_path}. Make sure it's a valid TSV file."
        all_pheno.append(error_message)
        all_scores.append(error_message)
        all_references.append(error_message)


    except Exception as e:
        error_message = f"An unexpected error occurred: {e}"
        all_pheno.append(error_message)
        all_scores.append(error_message)
        all_references.append(error_message)


clinvar_mono['Alpha_miss'] = all_pheno
clinvar_mono['Alpha_score'] = all_scores
clinvar_mono['aa_references'] = all_references

clinvar_mono.to_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_alphamiss.csv')

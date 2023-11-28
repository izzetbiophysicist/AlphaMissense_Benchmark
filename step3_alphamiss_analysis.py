#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 18:17:51 2023

@author: lucas
"""

####### Plot stuff
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ast
import os
import math
import re

os.chdir('/home/lucas/cataract/cataract_AlphaMiss/')
# Read the CSV file into a DataFrame
df = pd.read_csv("clinvar_alphamiss_annotation.csv")


##### Retrieve dates
df['date'] = df['Clinical significance (Last reviewed)'].str.extract(r'\((.*?)\)')
for i in range(len(df)):
    if pd.notna(df['date'][i]):
        print(i)
        
        df['date'][i] = df['date'][i].replace(" ", "")
        df['date'][i] = int(df['date'][i].split(',')[-1])
    else:
        df['date'][i] = 'no date'

df['Clinical significance (Last reviewed)'] = df['Clinical significance (Last reviewed)'].str.replace(r'\(.*?\)', '', regex=True)

# Define the segment to be removed
error_segment = 'Error: File not found at /home/lucas/cataract/cataract_AlphaMiss/database/P0DP91_output.tsv'

# Remove the segment from the entire column
df['Alpha_miss'] = df['Alpha_miss'].str.replace(error_segment, '')

###### Compress pathogenicity prediction for cases where multiple are found
alpha_list = []
for x in range(len(df['Alpha_miss'])):
    print(x)
    if(df['Alpha_miss'][x] != ''):
        tmp_list = ast.literal_eval(df['Alpha_miss'][x])
        if any(element == 'pathogenic' for element in tmp_list):
            result = 'Pathogenic'
        else:
            if any(element == 'benign' for element in tmp_list):
                result = 'Benign'
            if all(element == 'ambiguous' for element in tmp_list):
                result = 'Ambiguous'
                
        alpha_list.append(result)
    else:
        alpha_list.append('no prediction found')

df['Alpha_miss'] = alpha_list

### Retrieve review status
review_status = list(df['Review status'].unique())


def split_GO(input_string):
        
    # Use a regular expression to match the content inside brackets
    pattern = re.compile(r'\[([^]]+)\]')
    
    # Find all matches in the input string
    matches = pattern.findall(input_string)
    
    # Remove the content inside brackets from the input string
    result = re.sub(r'\[([^]]+)\]', '', input_string)
    
    # Split the result by semicolon and space
    result = [substring.strip() for substring in result.split(';')]
    
    return result

def accordance_by_GO(all_GO, filtered_df, column):
    ######
    filtered_df_noNAN = filtered_df[filtered_df[column].notna()]
    
    
    counts_series = pd.Series(all_GO).value_counts()
    
    # Convert the Series to a DataFrame
    result_df = counts_series.reset_index()
    result_df.columns = ['unique_GO', 'count']
    result_df['Accordance'] = 0
    result_df['Accordance_percent'] = 0
    result_df['unique_proteins'] = None
    
    for j in range(len(filtered_df_noNAN)):
        if filtered_df_noNAN.iloc[j]['Clinical significance (Last reviewed)'] == filtered_df_noNAN.iloc[j]['Alpha_miss']:
            current_GOs = split_GO(filtered_df_noNAN.iloc[j][column])
            for i in range(len(result_df)):
                if result_df.iloc[i]['unique_GO'] in current_GOs:
                    result_df.loc[i, 'Accordance'] += 1
    
    ### Count unique proteins
    for i in range(len(result_df)):
        print('Calculating for entry '+str(i)+' out of '+ str(len(result_df)))
        unique_proteins = []
        for j in range(len(filtered_df_noNAN)):
            if pd.isna(filtered_df_noNAN.iloc[j][column]) == False:
                current_GOs = split_GO(filtered_df_noNAN.iloc[j][column])
                
            if result_df.iloc[i]['unique_GO'] in current_GOs:
                unique_proteins.append(filtered_df_noNAN.iloc[j]['uniprot'])
        unique_proteins = len(set(unique_proteins))
        result_df.loc[i, 'unique_proteins'] = unique_proteins
            
            
    result_df['Accordance_percent'] = result_df['Accordance'] / result_df['count']                

    return result_df


def plot_prediction_by_review(df, review_st=None, date_cutoff=None, date_type='prev', approximate_clinvar=False, by_GO=False, plot_graph=True, classification=[]):
#### Filter dataframe
    filtered_df = df[df['Alpha_miss'] != 'no prediction found']
    
    if review_st != None:
        filtered_df = filtered_df[filtered_df['Review status'] == review_st]
    
    if date_cutoff != None and date_type=='prev':
        filtered_df = filtered_df[filtered_df['date'] != 'no date']
        filtered_df = filtered_df[filtered_df['date'] <= date_cutoff]
    
    if date_cutoff != None and date_type=='post':
        filtered_df = filtered_df[filtered_df['date'] != 'no date']
        filtered_df = filtered_df[filtered_df['date'] >= date_cutoff]
             
    if len(classification) != 0:
        filtered_df = filtered_df[filtered_df['Clinical significance (Last reviewed)'].isin(classification)]
    
    
    if approximate_clinvar:
        for i in range(len(filtered_df)):
            value = filtered_df.iloc[i]['Clinical significance (Last reviewed)']
            if value in ['Likely pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic']:
                filtered_df.iloc[i, filtered_df.columns.get_loc('Clinical significance (Last reviewed)')] = 'Pathogenic'
            
            if value in ['Benign/Likely benign', 'Likely benign', 'Benign']:
                filtered_df.iloc[i, filtered_df.columns.get_loc('Clinical significance (Last reviewed)')] = 'Benign'
    
    ### Calculate % accordance
    total_accord = 0
    for ground_truth, prediction in zip(filtered_df['Clinical significance (Last reviewed)'], filtered_df['Alpha_miss']):
        if ground_truth == prediction:
            total_accord += 1
        
    percent_accord = total_accord/len(filtered_df)
    
    #################
    #### Go counts
    #################
    
    
    all_GO_function = []
    all_GO_process = []
    all_GO_component = []
    

        
    for current_go in range(len(filtered_df)):
        if pd.isna(filtered_df.iloc[current_go]['GO_function']) == False:
            tmp = split_GO(filtered_df.iloc[current_go]['GO_function'])

            all_GO_function.extend(tmp)
            
            
        else:
            all_GO_function.append('None')
        
        if pd.isna(filtered_df.iloc[current_go]['GO_component']) == False:
            tmp = split_GO(filtered_df.iloc[current_go]['GO_component'])

            all_GO_component.extend(tmp)
        else:
            all_GO_component.append('None')
        
        if pd.isna(filtered_df.iloc[current_go]['GO_process']) == False:
            tmp = split_GO(filtered_df.iloc[current_go]['GO_process'])
            all_GO_process.extend(tmp)
            
        else:
            all_GO_process.append('None')
        
    #unique_function = set(all_GO_function)
    #unique_GO_function = list(unique_function)
    
    #unique_process = set(all_GO_process)
    #unique_GO_process = list(unique_process)
    
    
    #unique_component = set(all_GO_component)
    #unique_GO_component = list(unique_component)

    GOs = [all_GO_component, all_GO_function, all_GO_process]
    
    
    if plot_graph == True:
    
        # Set the order of the x-axis categories
        clinical_significance_order = [
            'not provided',
            'Uncertain significance/Uncertain risk allele',
            'Conflicting interpretations of pathogenicity',
            'Benign/Likely benign',
            'Likely benign',
            'Benign',
            'Uncertain significance',
            'Likely pathogenic',
            'Pathogenic/Likely pathogenic',
            'Pathogenic'
        ]
        
        # Set the color palette for each 'Alpha_miss' category
        palette = {'Benign': 'lightblue', 'Pathogenic': 'salmon', 'Ambiguous': 'lightgreen'}
        
        # Create a stacked bar plot with specified x-axis order
        plt.figure(figsize=(12, 6))
        
        # Count the occurrences
        counts = filtered_df.groupby(['Clinical significance (Last reviewed)', 'Alpha_miss']).size().unstack().fillna(0)
        
        # Normalize the counts
        normalized_counts = counts.div(counts.sum(axis=1), axis=0)
        
        # Plot the stacked bar plot
        normalized_counts.plot(kind='bar', stacked=True, color=[palette[col] for col in normalized_counts.columns], width=0.8)
        
    
        
        # Set plot labels and title
        plt.xlabel('Clinical Significance')
        plt.ylabel('Normalized Count')
        plt.title('Accordance '+str(cut))
        
        # Move the legend outside the plot area
        plt.legend(title='Alpha_miss', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Rotate x-axis labels for better visibility
        plt.xticks(rotation=45, ha='right')
        
        # Adjust layout for better display
        plt.tight_layout()
        
        # Show the plot
        plt.show()
    
    return percent_accord, len(filtered_df), filtered_df, GOs
 
############# Evaluate all non-VUS
classification= ['Benign/Likely benign',
'Likely benign',
'Benign',
'Likely pathogenic',
'Pathogenic/Likely pathogenic',
'Pathogenic']

accord_list = []   
sizes = []
for cut in range(2008, 2023):   
    result = plot_prediction_by_review(df, review_st=None, date_cutoff=cut, date_type='post', approximate_clinvar=True, plot_graph=True, classification=classification)
    accord_list.append(result[0])
    sizes.append(result[1])
    
############# Evaluate only pathogenic and benign
classification = ['Benign','Pathogenic']
accord_list_approximate = []
sizes_approximate = []
for cut in range(2008, 2023):   
    result = plot_prediction_by_review(df, review_st=None, date_cutoff=cut, date_type='post', approximate_clinvar=True, plot_graph=True, classification=classification)
    accord_list_approximate.append(result[0])
    sizes_approximate.append(result[1])

############# Evaluate only Non-determined
classification = ['not provided','Uncertain significance/Uncertain risk allele',
'Conflicting interpretations of pathogenicity',
'Uncertain significance']
accord_list_approximate = []
sizes_approximate = []
for cut in range(2008, 2023):   
    result = plot_prediction_by_review(df, review_st=None, date_cutoff=cut, date_type='post', approximate_clinvar=False, classification=classification)





############
###  Accordance by GO
############
classification = ['Benign','Pathogenic']
result = plot_prediction_by_review(df, review_st=None, date_cutoff=None, date_type='post', approximate_clinvar=True, plot_graph=False, classification=classification)
acc_Component = accordance_by_GO(result[3][0], result[2], column='GO_component')
acc_Function = accordance_by_GO(result[3][1], result[2], column='GO_function')
acc_Process = accordance_by_GO(result[3][2], result[2], column='GO_process')

acc_Component.to_csv('acc_Component.csv')
acc_Function.to_csv('acc_Function.csv')
acc_Process.to_csv('acc_Process.csv')


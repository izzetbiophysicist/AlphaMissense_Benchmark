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

def plot_prediction_by_review(df, review_st=None, date_cutoff=None, date_type='prev', strict=False):
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
    
    if strict == True:
        filtered_df = filtered_df[filtered_df['Clinical significance (Last reviewed)'].isin(['Pathogenic', 'Benign'])]
        
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
    plt.title('Stacked Bar Plot of Clinical Significance Colored by Alpha_miss Prediction')
    
    # Move the legend outside the plot area
    plt.legend(title='Alpha_miss', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Rotate x-axis labels for better visibility
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout for better display
    plt.tight_layout()
    
    # Show the plot
    plt.show()
    
    
plot_prediction_by_review(df, review_st=None, date_cutoff=2014, date_type='post', strict=True)
plot_prediction_by_review(df, review_st=None, date_cutoff=None, date_type='post', strict=True)

    
def plot_prediction_by_go():
    
    

def plot_prediction_by_pfam():
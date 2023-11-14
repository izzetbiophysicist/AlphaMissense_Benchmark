#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:24:08 2023

@author: lucas
"""

from unipressed import UniprotkbClient, UniparcClient
import time
from unipressed import IdMappingClient
import pandas as pd


def match_uniprot(clinvar, human_uniprot):
    
    all_matched = []
    pfam_list = []
    GO_function = []
    GO_component = []
    GO_process = []
    
    for g in clinvar['uniprot']:
        print('Finding uniprot ID for gene '+g)
        matches = []
        
        for h in range(len(human_uniprot)):
            if str(human_uniprot['Gene Names'][h]) != 'nan':
                if any(item in human_uniprot['Gene Names'][h].split() for item in g.split(';')):
                    ### Take the first match
                    matches = human_uniprot['Entry'][h]
                    break
                    
        #### Take the first one's annotation
        pfam_list.append(human_uniprot['Pfam'][h])
        GO_function.append(human_uniprot['Gene Ontology (molecular function)'][h])
        GO_component.append(human_uniprot['Gene Ontology (cellular component)'][h])
        GO_process.append(human_uniprot['Gene Ontology (biological process)'][h])                
        all_matched.append(matches)
        
    clinvar['Pfam'] = pfam_list
    clinvar['GO_function'] = GO_function
    clinvar['GO_component'] = GO_component
    clinvar['GO_process'] = GO_process
    
    return clinvar


clinvar = pd.read_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_alphamiss.csv')
human_uniprot = pd.read_csv('/home/lucas/cataract/cataract_AlphaMiss/uniprotkb_homo_sapiens_AND_model_organi_2023_11_01.tsv', sep='\t')

matched = match_uniprot(clinvar, human_uniprot)
matched.to_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_alphamiss_annotation.csv')



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


def retrieve_data(uniprot_result):
    
    uniprot_list = uniprot_result
    
    all_pfam_references = []
    all_interpro_references = []
    all_go_references = []
    all_reactome_references = []


    for uniprotID in uniprot_list:
        unikb = UniprotkbClient.fetch_one(uniprotID)
        cf=unikb.get('uniProtKBCrossReferences')
        pfam_references = []
        interpro_references = []
        go_references = []
        reactome_references = []
        
        print(uniprotID)
        
        for i in range(len(cf)):
            if cf[i].get('database') == 'Pfam':
                pfam_references.append(cf[i].get('id'))
    
            if cf[i].get('database') == 'GO':
                tmpgo=[]
                tmpgo=cf[i].get('properties')[0]
                go_references.append(tmpgo.get('value'))
                
                
            if cf[i].get('database') == 'InterPro':
                tmpip = []
                tmpip=cf[i].get('properties')[0]
                interpro_references.append(tmpip.get('value'))
    
            if cf[i].get('database') == 'Reactome':
                tmpreact=[]
                tmpreact=cf[i].get('properties')[0]
                reactome_references.append(tmpreact.get('value'))
                
        all_reactome_references.append(reactome_references)
        all_pfam_references.append(pfam_references)
        all_interpro_references.append(interpro_references)
        all_go_references.append(go_references)    
    
    gene_map = {
    #"Gene": uniprot_result['from'],
    "Uniprot": uniprot_list,
    "InterPro": all_interpro_references,
    "Pfam": all_pfam_references,
    "Go": all_go_references,
    "Reactome": all_reactome_references
    }
    
    gene_map_dataframe = pd.DataFrame(gene_map)
    
    return gene_map_dataframe


clinvar = pd.read_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_alphamiss.csv')
uniprot_ids = list(set(list(clinvar['uniprot'])))
gene_data = retrieve_data(uniprot_ids)

pfam = []
GO = []
for i in range(len(clinvar['uniprot'])):
    for j in range(len(gene_data['Uniprot'])):
        if clinvar['uniprot'][i] == gene_data['Uniprot'][j]:
                pfam.append(gene_data['Pfam'][j])
                GO.append(gene_data['Pfam'][j])

clinvar['Go'] = GO
clinvar['Pfam'] = pfam

clinvar.to_csv('/home/lucas/cataract/cataract_AlphaMiss/clinvar_alphamiss_annotation.csv')



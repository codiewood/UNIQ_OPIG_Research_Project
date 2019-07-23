#!/usr/bin/env python3
# coding: utf-8

#%%

import gzip, json, re, sys, os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

oas_file = str(sys.argv[1])
#oas_file = "Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz"
#oas_file = "Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz"

meta_line = True
sequence_data = []
for line in gzip.open(oas_file, 'rb'):
    if meta_line == True:
        metadata = json.loads(line)
        meta_line = False
        continue
    
    #Parse actual sequence data.
    basic_data = json.loads(line)
    sequence_data.append(basic_data)
    
    #IMGT-Numbered sequence.
    d = json.loads(basic_data['data'])
    sequence_data[-1]['data'] = d
    
#%%

def maximum_valued_key(dictionary):
    values = list(dictionary.values())
    keys = list(dictionary.keys())
    return keys[values.index(max(values))]

def sort_alphanumerically(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

    
#%%

unique_sequences = len(sequence_data)
regions = ['fwh1','cdrh1','fwh2','cdrh2','fwh3','cdrh3','fwh4']
amino_acid_list = ['A','G','I','L','P','V','F','W','Y','D','E','S','T','R','H','K','C','M','N','Q','Unused']

page = 1
for region_name in regions:
    def find_region_sequence(data):
        region_sequence = data['data'][region_name]
        return region_sequence

    def find_amino_acids(position):
        amino_acids = []
        for data in sequence_data:
            region_sequence = find_region_sequence(data)
            if position in region_sequence:
                amino_acids.append(region_sequence[position])
            else:
                amino_acids.append('Unused')
        return amino_acids
    
    def generate_positions():
        region_positions = set()
        for data in sequence_data:
            region_sequence = find_region_sequence(data)
            positions_used = list(region_sequence.keys())
            for position in positions_used:
                region_positions.add(position)
        region_position_list = sort_alphanumerically(list(region_positions))
        return region_position_list
    
    region_position_list = generate_positions()
    
    def amino_acid_percents(position):
        amino_acids = find_amino_acids(position)
        amino_acid_count = Counter(amino_acids)
        for amino_acid in amino_acid_count:
            frequency = (100*amino_acid_count[amino_acid])/unique_sequences
            amino_acid_count[amino_acid] = round(frequency,2)
        amino_acid_frequency = dict(amino_acid_count)
        return amino_acid_frequency
    
    raw_data = {}
    for amino_acid in amino_acid_list:
        raw_data[amino_acid] = []
        for position in region_position_list:
            amino_acid_frequency = amino_acid_percents(position)
            if amino_acid in amino_acid_frequency:
                raw_data[amino_acid].append(amino_acid_frequency[amino_acid])
            else:
                raw_data[amino_acid].append(0)
  
    df = pd.DataFrame(raw_data)
    colours = plt.cm.rainbow(np.linspace(0, 1, 21))
    colours[-1] = (0.95,0.95,0.95,1)
    
    species = str(metadata['Species']).replace('/','')
    
    sbc = df.plot.bar(stacked = True, color = colours,
                      title = str(metadata['Author'] + " " + metadata['Species']
                      + ' ' + region_name))
    sbc.set(xlabel="Position", ylabel="Frequency",
            xticklabels=region_position_list)
    sbc.legend(loc='center right', bbox_to_anchor=(-0.2, 0.5))
    sbc.get_figure().savefig('{}.pdf'.format(str(page) + "_" + species + "_"
                             + str(metadata['Size']) + "_" + region_name
                             + "_plot"), bbox_inches="tight")
    page += 1

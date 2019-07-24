#!/usr/bin/env python3
# coding: utf-8

#%%

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import oas_file_handling as oas

f = oas.oas_file(str(sys.argv[1]))
threshold = sys.argv[2] if len(sys.argv) > 2 else 1
#f = oas_file("Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz")
#f = oas_file("Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz")

amino_acid_list = f.amino_acids

page = 1
for region_name in f.regions:
    raw_data = {}
    for amino_acid in amino_acid_list:
        raw_data[amino_acid] = []
        for position in f.region_positions(region_name, freq, threshold):
            amino_acid_frequency = f.amino_acid_frequency(position)
            if amino_acid in amino_acid_frequency:
                raw_data[amino_acid].append(amino_acid_frequency[amino_acid])
            else:
                raw_data[amino_acid].append(0)
  
    df = pd.DataFrame(raw_data)
    colours = plt.cm.rainbow(np.linspace(0, 1, 21))
    colours[-1] = (0.95,0.95,0.95,1)
    
    species = str(f.metadata['Species']).replace('/','')
    
    sbc = df.plot.bar(stacked = True, color = colours,
                      title = str(f.metadata['Author'] + " " + f.metadata['Species']
                      + ' ' + region_name))
    sbc.set(xlabel="Position", ylabel="Frequency",
            xticklabels=f.region_positions(region_name, freq, threshold))
    sbc.legend(loc='center right', bbox_to_anchor=(-0.2, 0.5))
    sbc.get_figure().savefig('{}.pdf'.format(str(page) + "_" + species + "_"
                             + str(f.metadata['Size']) + "_" + region_name
                             + "_plot"), bbox_inches="tight")
    page += 1

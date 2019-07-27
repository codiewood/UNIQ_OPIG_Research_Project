#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import oas_file_handling as oas
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) > 2:
    length = int(sys.argv[2])
else:
    length = 15
f = oas.cdrh3_data(str(sys.argv[1]), length)

plot_data = {}
for amino_acid in f.amino_acids:
    plot_data[amino_acid] = []
    for position in f.region_positions('cdrh3'):
        if position in f.cdrh3_positions_used():
            amino_acid_frequency = f.cdrh3_amino_acid_frequency(position)
            if amino_acid in amino_acid_frequency:
                plot_data[amino_acid].append(amino_acid_frequency[amino_acid])
            else:
                plot_data[amino_acid].append(0)
     
plot_df = pd.DataFrame(plot_data)
colours = plt.cm.rainbow(np.linspace(0, 1, 21))
colours[-1] = (0.95,0.95,0.95,1)

sbc = plot_df.plot.bar(stacked = True, color = colours,
                       title = str(f.metadata['Author'] + " " + f.metadata['Species']
                  + " cdrh3 with lengths " + str(length)))
sbc.set(xlabel="Position", ylabel="Frequency",
        xticklabels=f.cdrh3_positions_used())
sbc.legend(loc='center right', bbox_to_anchor=(-0.2, 0.5))

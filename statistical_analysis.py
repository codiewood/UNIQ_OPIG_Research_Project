#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
from sklearn import metrics
from scipy.stats import chi2_contingency
from scipy.stats import chi2

f = oas.oas_file("Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz")
g = oas.oas_file("Collins_2015_IGHG_Mouse_sample_2_iglblastn_igblastn_IGHG.json.gz")
#g = oas.oas_file("Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz")
#f = oas.oas_file(str(sys.argv[1]))
#g = oas.oas_file(str(sys.argv[2]))

if len(sys.argv) > 3:
    threshold = int(sys.argv[3])
else:
    threshold = 1

def label_values(file,position):
    amino_acid_list = file.amino_acids
    label = []
    for amino_acid in amino_acid_list:
        amino_acid_count = file.amino_acid_occurences(position)
        if amino_acid in amino_acid_count:
            label.append(amino_acid_count[amino_acid])
        else:
            label.append(0)
    return label

labels_position1 = label_values(f,'1')
labels_position2 = label_values(f,'1')
lel = metrics.mutual_info_score(labels_position1, labels_position2)
lol = metrics.adjusted_mutual_info_score(labels_position1, labels_position2)
lul = metrics.normalized_mutual_info_score(labels_position1, labels_position2)
print(labels_position1)
print(labels_position2)
print(lel)
print(lol)
print(lul)
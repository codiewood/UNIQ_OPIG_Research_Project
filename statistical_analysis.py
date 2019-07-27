#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
import sys
from scipy.stats import chi2_contingency, chi2, spearmanr

#human_file = oas.oas_file("Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz")
#mouse_file = oas.oas_file("Collins_2015_IGHG_Mouse_sample_2_iglblastn_igblastn_IGHG.json.gz")
#mouse_file = oas.oas_file("Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz")

human_file = oas.oas_file(str(sys.argv[1]))
mouse_file = oas.oas_file(str(sys.argv[2]))

if len(sys.argv) > 3:
    threshold = int(sys.argv[3])
else:
    threshold = 1
    
if len(sys.argv) > 4:
    significance = int(sys.argv[4])
else:
    significance = 0.05

position_list = human_file.all_positions(threshold)

print("""Using data from {} and {} and a significance level of {}:
    """.format(human_file.file_name, mouse_file.file_name, significance))
for position in position_list:
    human_data = human_file.data_row(position)
    mouse_data = mouse_file.data_row(position)
    
    for i, amino_acid in reversed(list(enumerate(human_file.amino_acids))):
        if mouse_data[i] == 0 and human_data[i] == 0:
            del mouse_data[i]
            del human_data[i]
    
    table = []
    table.append(human_data)
    table.append(mouse_data)
    
    stat, p, dof, expected = chi2_contingency(table)
    
    prob = 1 - significance
    critical = chi2.ppf(prob, dof)
    print("""At position {}:
        """.format(position))
    print('Chi-Squared test:')
    print('Critical value = {}, Test statistic = {}'.format(critical, stat))
    print('Significance level = {}, p value = {}'.format(significance, p))
    if p <= significance:
    	print('We reject the null hypothesis that the human and mouse amino acid frequencies show no significant difference.')
    else:
    	print('We fail to reject the null hypothesis that the human and mouse amino acid frequencies show no significant difference.')
        
    coef, p = spearmanr(human_data, mouse_data)
    print('Spearmans Rank test:')
    print('Spearmans correlation coefficient: {}'.format(coef))
    if p <= significance:
    	print('We reject the null hypothesis that the human and mouse amino acid frequencies show no monotonic association.')
    else:
    	print('We fail to reject the null hypothesis that the human and mouse amino acid frequencies show no monotonic association.')
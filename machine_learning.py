#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt
import numpy as np

#%%
length = 15
trees = 10

human = oas.cdrh3_data('Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz', length)
mouse = oas.cdrh3_data('Collins_2015_IGHG_Mouse_sample_2_iglblastn_igblastn_IGHG.json.gz', length)
family = False

all_data = [human, mouse]

#%% Machine learning

if len(human.cdrh3_positions_used()) >= len(mouse.cdrh3_positions_used()):
    positions = human.cdrh3_positions_used()
else:
    positions = mouse.cdrh3_positions_used()

raw_data = {'Species': []}
for data in all_data:
    for k in range(data.number):
        if data.metadata['Species'] == 'human':
            raw_data['Species'].append(1)
        else:
            raw_data['Species'].append(0)
    for position in positions:
        for amino_acid in data.amino_acids:
            key = 'Position: '+ position +' AA: ' + amino_acid
            for cdrh3 in data.non_redundant:
                if key not in raw_data:
                    if position not in cdrh3:
                        if amino_acid == 'Unused':
                            raw_data[key] = [1]
                    elif cdrh3[position] == amino_acid:
                        raw_data[key] = [1]
                    else:
                        raw_data[key] = [0]
                else:
                    if position not in cdrh3:
                        if amino_acid == 'Unused':
                                raw_data[key].append(1)
                    elif cdrh3[position] == amino_acid:
                        raw_data[key].append(1)
                    else:
                        raw_data[key].append(0)

df = pd.DataFrame(raw_data)
attributes = df.iloc[:, 1:].values
labels = df.iloc[:, 0].values

training_attr, test_attr, training_labs, test_labs = train_test_split(attributes, labels, test_size=0.2, random_state=0)

RF = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
            max_depth=None, max_features='auto', max_leaf_nodes=None,
            min_impurity_decrease=0.0, min_impurity_split=None,
            min_samples_leaf=1, min_samples_split=2,
            min_weight_fraction_leaf=0.0, n_estimators=trees, n_jobs=1,
            oob_score=False, random_state=None, verbose=0,
            warm_start=False)

RF.fit(training_attr, training_labs)
predicted_labs = RF.predict(test_attr)
print("Accuracy:",metrics.accuracy_score(test_labs, predicted_labs))

position_importance = pd.Series(RF.feature_importances_)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(position_importance)


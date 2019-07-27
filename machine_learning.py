#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.model_selection import train_test_split
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics

#%%
length = 15
trees = 10

human = oas.cdrh3_data('Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz', length)
mouse = oas.cdrh3_data('Collins_2015_IGHG_Mouse_sample_2_iglblastn_igblastn_IGHG.json.gz', length)

all_data = [human, mouse]

#%% Machine learning

raw_data = {'Species': []}
for data in all_data:
    for k in range(data.number):
        if data.metadata['Species'] == 'human':
            raw_data['Species'].append(1)
        else:
            raw_data['Species'].append(0)
    for position in data.cdrh3_positions_used():
        key = 'Position '+ position
        for cdrh3 in data.sequences:
            if key not in raw_data:
                raw_data[key] = [cdrh3[position]]
            else:
                raw_data[key].append(cdrh3[position])

df = pd.DataFrame(raw_data)

for i, amino_acid in enumerate(human.amino_acids):
    df = df.replace(amino_acid,i)

attributes = df.iloc[:, 1:16].values
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

position_importance = pd.Series(RF.feature_importances_,index=human.cdrh3_positions_used())
print(position_importance)

bc = position_importance.plot.bar(color = colours, title = "Importance of positions")
bc.set(ylabel = 'Feature Importance Score', xlabel='Position')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:19:43 2019

@author: stats11
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

r = [0,1,2,3]
raw_data = {'1': [20, 5, 3, 2], '2': [5, 15, 5, 4],'3': [3, 5, 13, 17], '4': [2, 3, 7, 5], '5':[5, 7, 7, 7]}
df = pd.DataFrame(raw_data)
event_colours = plt.cm.rainbow(np.linspace(0, 1, 5))
sbc = df.plot.bar(stacked = True, color = event_colours)
sbc.set(xlabel="Position", ylabel="Frequency")
sbc.legend(loc='center right', bbox_to_anchor=(-0.2, 0.5))
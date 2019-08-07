#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
import sys

f = oas.oas_file(str(sys.argv[1]))

if len(sys.argv) > 2:
    family = str(sys.argv[2])
else:
    family = False

if len(sys.argv) > 3:
    threshold = int(sys.argv[3])
else:
    threshold = 1

print("""Using data from {}:
    """.format(f.file_name))
position_list = f.all_positions(family, threshold)
positions_used = []
for position1 in position_list:
    for position2 in position_list:
        if position2 not in positions_used:
            MI = f.mutual_information(position1, position2, family, normalized=True)
            if position1 == position2:
                print("""The self information of position {} has value {}.
                      """.format(position1,MI))
            else:
                print("""The mutual information of positions {} and {} has value {}.
                      """.format(position1,position2,MI))
    positions_used.append(position1)
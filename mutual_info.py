#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import oas_file_handling as oas
import sys

f = oas.oas_file(str(sys.argv[1]))

if len(sys.argv) > 2:
    threshold = int(sys.argv[2])
else:
    threshold = 1

print("""Using data from {}:
    """.format(f.file_name))
position_list = f.all_positions(threshold)
positions_used = []
for position1 in position_list:
    for position2 in position_list:
        if position2 not in positions_used:
            MI = f.mutual_information(position1, position2)
            if position1 == position2:
                print("""The self information of position {} has value {}.
                      """.format(position1,MI))
            else:
                print("""The mutual information of positions {} and {} has value {}.
                      """.format(position1,position2,MI))
    positions_used.append(position1)
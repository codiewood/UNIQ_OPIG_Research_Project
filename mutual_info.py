import sys, math
import oas_file_handling as oas

f = oas.oas_file("Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz")
#f = oas.oas_file("Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz")
#f = oas.oas_file(str(sys.argv[1]))
if len(sys.argv) > 2:
    threshold = int(sys.argv[2])
else:
    threshold = 1

def probability(position,amino_acid):
    amino_acid_frequency = f.amino_acid_frequency(position)
    if amino_acid in amino_acid_frequency:
        prob = amino_acid_frequency[amino_acid]/100
    else:
        prob = 0
    return prob
            
def joint_probability(position1,amino_acid1,position2,amino_acid2):
    both= 0
    for data in f.sequence_data:
        full_sequence = f.combined_sequence(data)
        if amino_acid2 == 'Unused' and position2 not in full_sequence:
            if amino_acid1 == 'Unused' and position1 not in full_sequence:
                both += 1
            elif position1 in full_sequence and full_sequence[position1] == amino_acid1:
                both += 1
        elif position2 in full_sequence and full_sequence[position2] == amino_acid2:
            if amino_acid1 == 'Unused' and position1 not in full_sequence:
                both += 1
            elif position1 in full_sequence and full_sequence[position1] == amino_acid1:
                both += 1
    joint_prob = both/f.unique_sequences
    return joint_prob

def mutual_information(position1,position2):
    MI = 0
    aa_position1 = set(f.find_amino_acids(position1))
    aa_position2 = set(f.find_amino_acids(position2))
    for amino_acid1 in aa_position1:
        for amino_acid2 in aa_position2:
            jp = joint_probability(position1,amino_acid1,position2,amino_acid2)
            p1 = probability(position1,amino_acid1)
            p2 = probability(position2,amino_acid2)
            if jp == 0:
                MI += 0
            else:
                MI += jp*math.log(jp/(p1*p2))
    MI_rounded = round(MI,4)
    return MI_rounded

#%%
print("""Using data from {}:
    """.format(f.file_name))
position_list = f.all_positions(threshold)
for position1 in position_list:
    for position2 in position_list:
        MI = mutual_information(position1,position2)
        if position1 == position2:
            print("""The self information of position {} has value {}.
                  """.format(position1,MI))
        else:
            print("""The mutual information of positions {} and {} has value {}.
                  """.format(position1,position2,MI))

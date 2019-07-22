import gzip, json, re, sys, math, os
from collections import Counter

#oas_file = str(sys.argv[1])
#oas_file = "Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz"
oas_file = "Corcoran_2016_heavy_mouse_IGHG_mouse_heavy_M2_igblastn_igblastn_IGHG.json.gz"

meta_line = True
sequence_data = []
for line in gzip.open(oas_file, 'rb'):
    if meta_line == True:
        metadata = json.loads(line)
        meta_line = False
        continue
    
    #Parse actual sequence data.
    basic_data = json.loads(line)
    sequence_data.append(basic_data)
    
    #IMGT-Numbered sequence.
    d = json.loads(basic_data['data'])
    sequence_data[-1]['data'] = d
    
def maximum_valued_key(dictionary):
    values = list(dictionary.values())
    keys = list(dictionary.keys())
    return keys[values.index(max(values))]

def sort_alphanumerically(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

redundant_sequences = metadata['Size']
regions = ['fwh1','cdrh1','fwh2','cdrh2','fwh3','cdrh3','fwh4']
amino_acid_list = ['A','G','I','L','P','V','F','W','Y','D','E','S','T','R','H','K','C','M','N','Q','Unused']

def combine_regions(data):
    combined_sequence = {}
    for region_name in regions:
        region_sequence = data['data'][region_name]
        combined_sequence.update(region_sequence)
    return combined_sequence 

def find_amino_acids(position):
    amino_acids = []
    for data in sequence_data:
        full_sequence = combine_regions(data)
        for i in range(int(data['redundancy'])):
            if position in full_sequence:
                amino_acids.append(full_sequence[position])
            else:
                amino_acids.append('Unused')
    return amino_acids

def generate_positions():
    all_positions = set()
    for data in sequence_data:
        full_sequence = combine_regions(data)
        positions_used = list(full_sequence.keys())
        for position in positions_used:
            all_positions.add(position)
    position_list = sort_alphanumerically(list(all_positions))
    return position_list 

position_list = generate_positions()

def amino_acid_percent(position):
    amino_acids = find_amino_acids(position)
    amino_acid_count = Counter(amino_acids)
    for amino_acid in amino_acid_count:
        frequency = (100*amino_acid_count[amino_acid])/redundant_sequences
        amino_acid_count[amino_acid] = round(frequency,2)
    amino_acid_frequency = dict(amino_acid_count)
    return amino_acid_frequency

def probability(position,amino_acid):
    amino_acid_frequency = amino_acid_percent(position)
    if amino_acid in amino_acid_frequency:
        prob = amino_acid_percent(position)[amino_acid]/100
    else:
        prob = 0
    return prob

def conditional_probability(X,x,Y,y): #prob that X=x given Y=y; X,Y positions and x,y amino acids
    condition = 0
    both = 0
    for data in sequence_data:
        full_sequence = combine_regions(data)
        if Y in full_sequence and full_sequence[Y] == y:
            condition += 1
            if X in full_sequence and full_sequence[X] == x:
                both += 1
    if condition == 0:
        cond_prob = 0
    else:
        cond_prob = both/condition
    return cond_prob
            
def joint_probability(X,x,Y,y):
    joint_prob = probability(Y,y)*conditional_probability(X,x,Y,y)
    return joint_prob

def mutual_information(X,Y):
    MI = 0
    for amino_acid1 in amino_acid_list:
        for amino_acid2 in amino_acid_list:
            jp = joint_probability(X,amino_acid1,Y,amino_acid2)
            px = probability(X,amino_acid1)
            py = probability(Y,amino_acid2)
            if px*py == 0 or jp == 0:
                MI += 0
            else:
                MI += jp*math.log(jp/(px*py))
    return MI        

print("""Using data from {}:
              """.format(oas_file))
for position1 in position_list:
    for position2 in position_list:
        MI = mutual_information(position1,position2)
        if position1 == position2:
            print("""The self information of position {} has value {}.
              """.format(position1,MI))
        else:
            print("""The mutual information of positions {} and {} has value {}.
              """.format(position1,position2,MI))

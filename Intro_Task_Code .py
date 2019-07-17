#!/usr/bin/env python
# coding: utf-8

# In[2]:


import gzip
import json

oas_file = "Vander_Heiden_2017_Heavy_HD09_IGHG_HD09_Unsorted_Bcells_age31_healthy_iglblastn_igblastn_IGHG.json.gz"

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


# In[3]:


unique_sequences = len(sequence_data)
print("There are {} unique sequences in this data.".format(unique_sequences))


# In[4]:


redundant_sequences = 0
for data in sequence_data:
    duplicates = data['redundancy']
    redundant_sequences += duplicates

print("There are {} redundant sequences in this data.".format(redundant_sequences))    


# In[5]:


max_errors = 0
correct_sequences = 0
for data in sequence_data:
    if data['num_errors'] == "0":
        correct_sequences += 1
    elif int(data['num_errors']) > max_errors:
        max_errors = int(data['num_errors'])
erroneous_sequences = unique_sequences - correct_sequences

print("There are {} sequences thought to contain errors in this data.".format(erroneous_sequences))
print("{} sequences are thought to contain 0 errors.".format(correct_sequences))

for i in range(1, max_errors + 1):
    i_errors = 0
    for data in sequence_data:
        if int(data['num_errors']) == i:
            i_errors += 1
    print("{} sequences are thought to contain {} errors.".format(i_errors,i))


# In[6]:


cdrh3_sequences = []
for data in sequence_data:
    cdrh3 = data['data']['cdrh3']
    cdrh3_sequences.append(cdrh3)

cdh3_sequence_counter = {}
for cdrh3 in cdrh3_sequences:
    key = json.dumps(cdrh3, sort_keys=True)
    if key in cdh3_sequence_counter:
        cdh3_sequence_counter[key] += 1
    else:
        cdh3_sequence_counter[key] = 1

def maximum_valued_key(dictionary):
    values = list(dictionary.values())
    keys = list(dictionary.keys())
    return keys[values.index(max(values))]

most_common = maximum_valued_key(cdh3_sequence_counter)
occurences = cdh3_sequence_counter[most_common]

print("The CDR-H3 sequence which was most common was {}, which appeared {} times.".format(most_common,occurences))
    


# In[7]:


import re

def combine_regions(data):
    regions = ['fwh1','cdrh1','fwh2','cdrh2','fwh3','cdrh3','fwh4']
    combined_sequence = {}
    for region_name in regions:
        region_sequence = data['data'][region_name]
        combined_sequence.update(region_sequence)
    return combined_sequence #dictionary of the whole sequence, not split into regions

def find_amino_acids(position):
    amino_acids = []
    for data in sequence_data:
        full_sequence = combine_regions(data)
        if position in full_sequence:
            amino_acids.append(full_sequence[position])
    return amino_acids #list of amino acids at that position across all data (inc repeats)

def sort_alphanumerically(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def generate_positions():
    all_positions = set()
    for data in sequence_data:
        full_sequence = combine_regions(data)
        positions_used = list(full_sequence.keys())
        for position in positions_used:
            all_positions.add(position)
    position_list = sort_alphanumerically(list(all_positions))
    return position_list #list of all possible positions in order of occurence

def percentage_usage(amino_acid):
    percent = (100*amino_acid_count[amino_acid])/unique_sequences
    rounded_percent = round(percent,2)
    return rounded_percent 
    
consensus_sequence = {}
position_list = generate_positions()
for position in position_list:
    amino_acids = find_amino_acids(position)
    amino_acid_count = Counter(amino_acids)
    for amino_acid in amino_acid_count:
        amino_acid_count[amino_acid] = percentage_usage(amino_acid)
    amino_acid_usage = dict(amino_acid_count)
    
    print("""For position {}, we see that the amino acids seen and their respective percentage usages are:
    {}.
    """.format(position,amino_acid_usage))

    consensus_sequence[position] = maximum_valued_key(amino_acid_usage)

print("""The consensus sequence for this data is:
{}""".format(consensus_sequence))


# In[8]:


def find_amino_acids_repeats(position):
    amino_acids = []
    for data in sequence_data:
        full_sequence = combine_regions(data)
        if position in full_sequence:
            for i in range(int(data['redundancy'])):
                amino_acids.append(full_sequence[position])
    return amino_acids #list of amino acids at that position across all data (inc repeats)

def percentage_usage_repeats(amino_acid):
    percent = (100*amino_acid_count[amino_acid])/redundant_sequences
    rounded_percent = round(percent,2)
    return rounded_percent 
    
consensus_sequence = {}
position_list = generate_positions()
for position in position_list:
    amino_acids = find_amino_acids_repeats(position)
    amino_acid_count = Counter(amino_acids)
    for amino_acid in amino_acid_count:
        amino_acid_count[amino_acid] = percentage_usage_repeats(amino_acid)
    amino_acid_usage = dict(amino_acid_count)
    
    print("""For position {}, we see that the amino acids seen and their respective percentage usages, including sequence repeats, are:
    {}.
    """.format(position,amino_acid_usage))

    consensus_sequence[position] = maximum_valued_key(amino_acid_usage)

print("""The consensus sequence for this data, including sequence repeats, is:
{}""".format(consensus_sequence))


# In[ ]:





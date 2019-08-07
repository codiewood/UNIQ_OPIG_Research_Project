#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip, json, re
from collections import Counter
from sklearn import metrics

def maximum_valued_key(dictionary):
    values = list(dictionary.values())
    keys = list(dictionary.keys())
    return keys[values.index(max(values))]

def sort_alphanumerically(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)

class oas_file():
    def __init__(self, src): 
        meta_line = True
        seq_data = []
        for line in gzip.open(src,'rb'):
            if meta_line == True:
                    self.metadata = json.loads(line) 
                    meta_line=False
                    continue
            #Parse actual sequence data.
            basic_data = json.loads(line)
            seq_data.append(basic_data)
            
            #IMGT-Numbered sequence.
            d = json.loads(basic_data['data'])
            seq_data[-1]['data'] = d
        self.sequence_data = seq_data
        if self.metadata['Chain'] == 'Heavy':
            self.regions = ['fwh1','cdrh1','fwh2','cdrh2','fwh3','cdrh3','fwh4']
        elif self.metadata['Chain'] == 'Light':
            self.regions = ['fwl1','cdrl1','fwl2','cdrl2','fwl3','cdrl3','fwl4']
        self.amino_acids = ['A','G','I','L','P','V','F','W','Y','D','E','S','T','R','H','K','C','M','N','Q','Unused']
        self.file_name = src
        gene_annotations = set()
        for seq in self.sequence_data:
            gene_annotations.add(seq['v'])
        self.gene_annotations = list(gene_annotations)
        fams = set()
        for text in self.gene_annotations:
            fam = text.partition("-")[0]
            if 'S' in fam:
                fam = fam.partition("S")[0]
            fams.add(fam)
        self.families = list(fams)
    
    def error_counts(self):
        max_errors = 0
        correct_sequences = 0
        for data in self.sequence_data:
            if data['num_errors'] == "0":
                correct_sequences += 1
            elif int(data['num_errors']) > max_errors:
                max_errors = int(data['num_errors'])
        erroneous_sequences = self.unique_sequences - correct_sequences
        
        print("There are {} sequences thought to contain errors in this data."
              .format(erroneous_sequences))
        print("{} sequences are thought to contain 0 errors."
              .format(correct_sequences))
        
        for i in range(1, max_errors + 1):
            i_errors = 0
            for data in self.sequence_data:
                if int(data['num_errors']) == i:
                    i_errors += 1
            print("{} sequences are thought to contain {} errors."
                  .format(i_errors,i))            
    
    def family_sequences(self, family=False):
        family_sequences = []
        if family == False:
            return self.sequence_data
        else:
            for data in self.sequence_data:
                data_family = data['v'].partition("-")[0]
                if 'S' in data_family:
                    data_family = data_family.partition("S")[0]
                if data_family == family:
                    family_sequences.append(data)
        return family_sequences
    
    def unique_sequences(self, family=False):
        family_sequences = self.family_sequences(family)
        num = len(family_sequences)
        return num
    
    def region_sequences(self,region_name, family=False):
        region_sequences = []
        for data in self.family_sequences(family):
            region = data['data'][region_name]
            region_sequences.append(region)
        return region_sequences
    
    def combined_sequence(self, data):
        combined_sequence = {}
        for region_name in self.regions:
            region_sequence = data['data'][region_name]
            combined_sequence.update(region_sequence)
        return combined_sequence
    
    def most_common_sequence(self, region_name, family=False):
        reg_sequence_counter = {}
        for seq in self.region_sequences(region_name,family):
            key = json.dumps(seq, sort_keys=True)
            if key in reg_sequence_counter:
                reg_sequence_counter[key] += 1
            else:
                reg_sequence_counter[key] = 1
        most_common = maximum_valued_key(reg_sequence_counter)
        occurences = reg_sequence_counter[most_common]
        if family == False:
            print("The {} sequence which was most common was {}, which appeared {} times."
                  .format(region_name,most_common,occurences))
        else:
            print("The {} sequence in family {} which was most common was {}, which appeared {} times."
                  .format(region_name,family,most_common,occurences))
        
    def region_positions(self, region_name, family=False, threshold=1):
        region_positions = set()
        for data in self.region_sequences(region_name,family):
                positions_used = list(data.keys())
                for position in positions_used:
                    region_positions.add(position)
        region_position_list = sort_alphanumerically(list(region_positions))
        if threshold > 1:
            freq_positions = []
            for position in region_position_list:
                position_use_count = 0
                for data in self.region_sequences(region_name,family):
                    if position in data:
                        position_use_count += 1
                if position_use_count >= threshold:
                    freq_positions.append(position)
            return freq_positions
        else:
            return region_position_list
    
    def all_positions(self, family=False, threshold=1):
        all_positions = set()
        for region_name in self.regions:
            for position in self.region_positions(region_name, family, threshold):
                all_positions.add(position)
        all_positions_list = sort_alphanumerically(list(all_positions))
        return all_positions_list
    
    def position_use_count(self, position, region_name, family=False):
        position_use_count = 0
        for data in self.region_sequences(region_name,family):
            if position in data:
                position_use_count += 1
        return position_use_count
    
    def find_amino_acids(self, position, family=False):
        amino_acids = []
        for data in self.family_sequences(family):
            if position in self.combined_sequence(data):
                amino_acids.append(self.combined_sequence(data)[position])
            else:
                amino_acids.append('Unused')
        return amino_acids
        
    def amino_acid_occurences(self, position, family=False):
        amino_acids = self.find_amino_acids(position, family)
        amino_acid_count = dict(Counter(amino_acids))
        return amino_acid_count

    def amino_acid_frequency(self, position, family=False):
        amino_acid_count = self.amino_acid_occurences(position,family)
        for amino_acid in amino_acid_count:
            frequency = (100*amino_acid_count[amino_acid])/(self.unique_sequences(family))
            amino_acid_count[amino_acid] = round(frequency,2)
        return amino_acid_count
    
    def print_frequencies(self, family=False, threshold=1):
        freq_dict = {}
        for position in self.all_positions(family, threshold):
            freq_dict[position] = self.amino_acid_frequency(position,family)
        return freq_dict

    def consensus_sequence(self, family=False, threshold=1):
        consensus_sequence = {}
        for position in self.all_positions(family,threshold):
            consensus_sequence[position] = maximum_valued_key(self.amino_acid_frequency(position,family))
        print("""The consensus sequence for this data is:
            {}""".format(consensus_sequence))
            
    def data_row(self,position,family=False):
        amino_acid_list = self.amino_acids
        row = []
        for amino_acid in amino_acid_list:
            amino_acid_count = self.amino_acid_occurences(position,family)
            if amino_acid in amino_acid_count:
                row.append(amino_acid_count[amino_acid])
            else:
                row.append(0)
        return row
    
    def mutual_information(self, position1, position2, family=False, normalized=False, adjusted=False):
        labels_position1 = self.data_row(position1,family)
        labels_position2 = self.data_row(position2,family)
        if normalized == True:
            MI = metrics.normalized_mutual_info_score(labels_position1, labels_position2)
        elif adjusted == True:
            MI = metrics.adjusted_mutual_info_score(labels_position1, labels_position2)
        else:
            MI = metrics.mutual_info_score(labels_position1, labels_position2)
        return MI
    

class cdrh3_data(oas_file):
    def __init__(self, src, length, family=False):
        super(cdrh3_data, self).__init__(src)
        cdrh3_sequences = self.region_sequences('cdrh3',family)
        cdrh3_same_length = []
        for cdrh3 in cdrh3_sequences:
            if len(cdrh3) == length:
                cdrh3_same_length.append(cdrh3)
        self.sequences = cdrh3_same_length
        self.family = family
        self.length = length
        self.redundant = len(self.sequences)
        self.non_redundant = [dict(t) for t in {tuple(seq.items()) for seq in self.sequences}]
        self.number = len(self.non_redundant)        
    
    def find_cdrh3_amino_acids(self,position):
        amino_acids = []
        for cdrh3 in self.non_redundant:
            if position in cdrh3:
                amino_acids.append(cdrh3[position])
            else:
                amino_acids.append('Unused')
        return amino_acids
        
    def cdrh3_amino_acid_occurences(self, position):
        amino_acids = self.find_cdrh3_amino_acids(position)
        amino_acid_count = dict(Counter(amino_acids))
        return amino_acid_count
    
    def cdrh3_amino_acid_frequency(self, position):
        amino_acid_count = self.cdrh3_amino_acid_occurences(position)
        for amino_acid in amino_acid_count:
            frequency = (100*amino_acid_count[amino_acid])/self.number
            amino_acid_count[amino_acid] = round(frequency,2)
        return amino_acid_count
    
    def cdrh3_position_use_count(self, position):
        position_use_count = 0
        for cdrh3 in self.non_redundant:
            if position in cdrh3:
                position_use_count += 1
        return position_use_count
    
    def cdrh3_positions_used(self):
        positions_used = []
        for position in self.region_positions('cdrh3', self.family):
                if self.cdrh3_position_use_count(position) >= 1:
                    positions_used.append(position)
        return positions_used
    
    
            
            
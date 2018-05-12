#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 13:15:37 2018

@author: martin
"""
import time
from collections import defaultdict


start = time.clock()


include_num_replicates = 7

class Protein(object):
    
    def __init__(self, accession, condition, replicate, quantity):
        self.accession = accession
        self.condition = condition
        self.replicate = replicate
        self.quantity = quantity

pdict = defaultdict(list)

conds = []

# def parse():
with open('EF.MM.Body_donors.OA.180404.txt', 'r') as report:
    next(report)
    for line in report:
        line = line.strip()
        line = line.split('\t')
        condition = line[0]
        if condition not in conds:
            conds.append(condition)
        replicate = int(line[2])
        #if replicate not in reps: reps.add(replicate)
        
        accession = line[3].split(';')[0]
        
        quantity = line[4]
        quantity = quantity.replace('"','')
        quantity = quantity.replace(",", ".")
        
        

        if replicate <= include_num_replicates:
            prot = Protein(accession, condition, replicate, quantity)
            pdict[prot.accession].append(prot)


# For every accession number as key, initiate list of length (num conditions * num replicates)
pmatrix = {p:[0]*len(conds)*include_num_replicates for p in pdict.keys()}



# def populate_matrix():
for acc, plist in pdict.items():
    for protein in plist:
        base_value = conds.index(protein.condition) * include_num_replicates
        protein_index = base_value + (protein.replicate - 1)

        pmatrix[protein.accession][protein_index] = protein.quantity



with open('OA.donors.body.full.wide.csv', 'w') as fw:
    header = ','.join(["%s.%s"%(x,y) for x in conds for y in range(1,include_num_replicates+1)])
    fw.write('Accession')
    fw.write(',')
    fw.write(header)
    fw.write('\n')
    for key, value in pmatrix.items():
        value = ','.join([str(x) for x in value])
        fw.write(key)
        fw.write(',')
        fw.write(value)
        fw.write('\n')


print(time.clock() - start)

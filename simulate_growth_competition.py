#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:55:28 2018

Simulate a CRISPR screen, starting from a sgRNA library. 
Transduction, selection and sequencign are simulated, leading to a final output of sgRNA counts.

Parameters:
    --cov_virus    transduction coverage
    --cov_cells    cell culture coverage
    --cov_pcr      sequencing coverage
    --dupl_time    cell duplication time
    --freq_negfc   frequency of essential (lethal)
    --freq_posfc   frequency of growth supressing
    --n_repl_lib_pcr number of sequencing replicates for the library
    --n_repl_sel   number of selection phase replicates
    --n_repl_pcr   number of sequencing replicates for selected timepoints
    --n_splittings number of cell splittings during selection phase
    --n_sgrnas     number of sgRNAs in total
    --n_sgrnas_per_gene  number of sgRNAs per gene
    --lib_breadth  library breadth

@author: Katharina Imkeller
"""

import numpy as np
import argparse as argparse


### PARSE arguments

parser = argparse.ArgumentParser()
parser.add_argument('--cov_virus', 
                    default = 400, type = int, 
                    help='Coverage during viral transduction.')
parser.add_argument('--cov_cells', 
                    default = 400, type = int, 
                    help='Coverage during cell culture selection.')
parser.add_argument('--cov_pcr', 
                    default = 400, type = int, 
                    help='Coverage for PCR amplification.')
parser.add_argument('--dupl_time', 
                    default = 30, type = int, 
                    help='Cell duplication time in hours.')
parser.add_argument('--freq_negfc', 
                    default = 0.1, type = float, 
                    help='Frequency of sgRNAs with negative fitness effect.')
parser.add_argument('--freq_posfc', 
                    default = 0.01, type = float, 
                    help='Frequency of sgRNAs with positive fitness effect.')
parser.add_argument('--n_repl_lib_pcr', 
                    default = 2, type = int, 
                    help='Number of replicates selection.')
parser.add_argument('--n_repl_sel', 
                    default = 10, type = int, 
                    help='Number of replicates selection.')
parser.add_argument('--n_repl_pcr', 
                    default = 3, type = int, 
                    help='Number of replicates final PCR.')
parser.add_argument('--n_splittings', 
                    default = 7, type = int, 
                    help='Number of cell splittings.')
parser.add_argument('--n_sgrnas', 
                    default = 50000, type = int, 
                    help='Number of cell splittings.')
parser.add_argument('--n_sgrnas_per_gene', 
                    default = 4, type = int, 
                    help='Number of cell splittings.')
parser.add_argument('--lib_breadth', 
                    default = 1.1, type = float, 
                    help='Breadth of a the initial sgRNA library. \
                    Fix it to specific choices, because there is only a limited number of meaningfull choices. \
                    Can be between 0.5 and 1.5')
args = parser.parse_args()




### BUILD LIBRARY

# lognormal library ditribution
# array of counts per sgRNA
library_counts = np.random.lognormal(5, args.lib_breadth, size=args.n_sgrnas).astype(int)

# attribute a growth rate to every sgRNA
# some of the id counts are essential, some are tumorsupressors
p_essential = args.freq_negfc + args.freq_posfc
# select which ones will be affected
selected = np.random.choice(args.n_sgrnas, size=int(args.n_sgrnas*p_essential))
# this has to be done in two steps, so that one sgrna does not get picked for pos and neg
# one part will have neg FC (essentials)
# get indices
essentials_positions = selected[0:int(args.n_sgrnas*args.freq_negfc)]
# one part will have pos FC
growthsuppressors_positions = selected[int(args.n_sgrnas*args.freq_negfc):]

# calculate the growth rate from the duplication time
# one cell splitting lasts 3 days - 72 h
growth_rate = np.log(2**(72./args.dupl_time))

# build a growth rate array for every sgRNA
growthrate_array = np.empty((1, args.n_sgrnas))[0]
class_list = []
gene_list = []
# count sgRNAs to attribute gene ID
essential_counter = 0
growths_counter = 0
neutral_counter = 0

for index in np.arange(args.n_sgrnas):
    if index in essentials_positions:
        # growth rates are affected 0-100%
        growthrate_array[index] = growth_rate - growth_rate * np.random.choice(np.arange(0,0.21,0.01))
        class_list.append('essential')
        gene_list.append(str(int(essential_counter / args.n_sgrnas_per_gene)))
        essential_counter = essential_counter + 1
    elif index in growthsuppressors_positions:
        growthrate_array[index] = growth_rate + growth_rate * np.random.choice(np.arange(0,0.21,0.01))
        class_list.append('growthsupressing')
        gene_list.append(str(int(growths_counter / args.n_sgrnas_per_gene)))
        growths_counter = growths_counter + 1
    else: 
        growthrate_array[index] = growth_rate
        class_list.append('neutral')
        gene_list.append(str(int(neutral_counter /  args.n_sgrnas_per_gene)))
        neutral_counter = neutral_counter + 1



### FUNCTIONS
        
def sequencing_result(count_array):
    # goal: generate sequencign results for a given sgRNA pool
    # comment: PCR amplification is assumed to be linear, only pick-up likelihood is considered
    # calculate probability
    results = np.empty((1, args.n_sgrnas))[0]
    prob_pickup = args.n_sgrnas * args.cov_pcr / float(sum(count_array))
    if (prob_pickup < 1):
        for i in range(len(count_array)):
            new_count = sum(np.random.binomial(1, p=prob_pickup, size=int(count_array[i])))
            results[i] = new_count
        return results
    else:
        count_array = grow(count_array)
        return(sequencing_result(count_array))

def transduce(count_array):
    # modify the count array during transduction
    prob_pickup = args.n_sgrnas * args.cov_virus / float(sum(count_array))
    if (prob_pickup < 1):
        for i in range(len(count_array)):
            new_count = sum(np.random.binomial(1, p=prob_pickup, size=int(count_array[i])))
            count_array[i] = new_count
    else:
        pass

def split(count_array):
    # modify the count array during transduction
    prob_pickup = args.n_sgrnas * args.cov_cells / float(sum(count_array))
    if (prob_pickup < 1):
        for i in range(len(count_array)):
                new_count = sum(np.random.binomial(1, p=prob_pickup, size=int(count_array[i])))
                count_array[i] = new_count
    else:
        count_array = grow(count_array)
        split(count_array)

def grow(count_array):
    count_array = np.multiply(count_array, np.exp(growthrate_array)).astype(int)
    return(count_array)


# sequence the library n_repl_lib_pcr times
library_sequencing_counts = sequencing_result(library_counts)
for i in np.arange(args.n_repl_lib_pcr-1):
    library_sequencing_counts = np.vstack((library_sequencing_counts,
                                               sequencing_result(library_counts)))


# transduce the library ONCE
transduce(library_counts)
# grow ONCE
library_counts = grow(library_counts)
# sample timepoint T0, sequencing
T0_sequencing_counts = sequencing_result(library_counts)


### SELECTION

# list holding replicate_sequencings
replicate_results = []
# number of replicates
for simulation in np.arange(args.n_repl_sel):
    # generate new copy of the array for every replicate
    replicate_count = np.copy(library_counts)
    
    # number of cell culturing passages
    for passage in np.arange(args.n_splittings):
        # split
        split(replicate_count)
        # grow
        # why does it not automatically change the array?
        replicate_count = grow(replicate_count)
    
    # sequence
    replicate_sequencing_counts = sequencing_result(replicate_count)
    for i in np.arange(args.n_repl_pcr):
        replicate_sequencing_counts = np.vstack((replicate_sequencing_counts,
                                                sequencing_result(replicate_count)))
    
    replicate_results.append(replicate_sequencing_counts)
    

### writing OUTPUT
outfile_res = open('results_' +
               str(args.cov_cells) + "_" + 
               str(args.cov_pcr) + "_" + 
               str(args.cov_virus) + "_" + 
               str(args.dupl_time) + "_" + 
               str(args.freq_negfc) + "_" + 
               str(args.freq_posfc) + "_" + 
               str(args.lib_breadth) + "_" + 
               str(args.n_repl_lib_pcr) + "_" + 
               str(args.n_repl_sel) + "_" + 
               str(args.n_repl_pcr) + "_" + 
               str(args.n_sgrnas) + "_" + 
               str(args.n_splittings) + '.tsv', 'w')

outfile_res.write(str(args) + "\n" + "sgrna_id\t" +
              "growth_rate\t" + "class\t" + "gene_id" +
              "".join(["\tlibrary" + str(n) for n in np.arange(args.n_repl_lib_pcr)]) +
              "\tT0" + 
              "".join(["".join(["\tR" + str(n) + "_" + str(i) for i in np.arange(args.n_repl_pcr)])
                      for n in np.arange(args.n_repl_sel)]) +
              "\n")

for index in np.arange(args.n_sgrnas):
    outfile_res.write(str(index) + "\t")
    outfile_res.write(str(growthrate_array[index]) + "\t")
    outfile_res.write(str(class_list[index]) + "\t")
    outfile_res.write(str(gene_list[index]) + "\t")
    # write all library sequencing results
    for i in np.arange(args.n_repl_lib_pcr):
        outfile_res.write(str(library_sequencing_counts[i][index]) + "\t")
    outfile_res.write(str(T0_sequencing_counts[index]))
    for repl in np.arange(args.n_repl_sel):
        # write all replicate sequencing results
        for i in np.arange(args.n_repl_pcr):
            outfile_res.write("\t" + str(replicate_results[repl][i][index]))
    outfile_res.write("\n")

outfile_res.close()
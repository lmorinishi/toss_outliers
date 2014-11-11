#!/bin/python
from __future__ import division
import cPickle as pickle
import sys
import os
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import pprint

#Functions
'''
def calculate_fitness(barcode, counts, timepoints, labels=None):
    if len(counts) == 0:
        return None
    log_counts = np.log2(np.array(counts))
    regression = np.polyfit(timepoints, log_counts, deg=1)
    slope, intercept = regression[0:2]

    return slope

def calculate_fitness_paper(barcode, counts, wt_counts, timepoints, labels=None):
    if len(counts) == 0:
        return None
    counts_over_wildtype = np.array(counts) / np.array(wt_counts)
    log_counts = np.log2(counts_over_wildtype)
    regression = np.polyfit(timepoints, log_counts, deg=1)
    slope, intercept = regression[0:2]

    plt.scatter(timepoints, log_counts)
    fitted_y = slope*timepoints + intercept
    plt.plot(timepoints, fitted_y)
    plt.xlabel('time (s)')
    plt.ylabel('log2 read count')
    plt.show()
    return slope
'''

def barcodes_to_codons_fitness(barcode_fitness_scores, allele_dict):
    codon_fitness_dict = {}

    for barcode, mutant_tuple in allele_dict.items():
        pos, codon = mutant_tuple
        codon_fitness_dict.setdefault(pos, {})
        codon_fitness_dict[pos].setdefault(codon, {})
        if barcode in barcode_fitness_scores:
            codon_fitness_dict[pos][codon][barcode] = fitness

    return codon_fitness_dict

def rna_to_nt(rna):
    return rna.replace('U','T').replace('u','t')

def nt_to_rna(nt):
    return nt.replace('T','U').replace('t','u')

def filtered_mean(data):
    data_array = np.array(data)
    return np.mean(data_array)

def codons_to_aa_fitness(codon_fitness_dict, translate_dict):
    aa_fitness_dict = {}

    for pos in sorted(codon_fitness_dict.keys()):
        aa_fitness.setdefault(pos, {})
        for codon, aa in translate_dict.items():
            aa_fitness[pos].setdefault(aa, {})
            rna_codon = nt_to_rna(codon)

            if rna_codon in codon_fitness_dict[pos]:
                barcode_fitness_values = [x[1] for x in codon_fitness_dict[pos][rna_codon].items()]
                codon_fitness = filtered_mean(barcode_fitness_scores)
                if len(codon_fitness) == 0:
                    aa_fitness[pos][aa][codon] = None
                else:
                    aa_fitness[pos][aa][codon] = codon_fitness
            else:
                aa_fitness[pos][aa][codon] = None

    return aa_fitness_dict

def make_aa_fitness_matrix(aa_fitness_scores, amino_to_number_dict):
    fitness_matrix = np.empty(shape=(len(aa_fitness_scores),len(amino_to_number_dict)))

    i = 0
    for pos, aa_fitness in sorted(aa_fitness_scores.items()):
        for aa in aa_fitness:
            aa_num = amino_to_number_dict[aa]
            fitness_values = np.array([v for k,v in aa_fitness.items() if v is not None])
            if len(fitness_values) > 0:
                fitness_matrix[i][aa_num] = fitness_values
            else:
                fitness_matrix[i][aa_num] = None
        i+=1

    return fitness_matrix


#Get files from command line
(aa_fitness_pickle, amino_to_number_pickle) = sys.argv[1:2]

#Read in files
aa_fitness_dict = pickle.load(open(aa_fitness_pickle, "rb"))
amino_to_number_dict = pickle.load(open(amino_to_number_pickle, "rb"))

aa_fitness_matrix = make_aa_fitness_matrix(aa_fitness_dict)
#codon_to_number_dict = get_codon_to_number_dict(translate_dict, amino_to_number_dict)
#Plot heatmaps of barcode counts for each amino acid and position
column_labels = [x for x,y in sorted(amino_to_number_dict.items(), key=lambda x: x[1])]
row_labels = sorted(aa_fitness_scores.keys())
ax1 = plt.subplot(211) #this plot will use a linear color scale
heatmap = plt.imshow(aa_fitness_matrix, cmap=plt.cm.YlGnBu, interpolation="nearest")
ax1.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
#ax1.set_yticks(np.arange(barcode_count_aa_heatmap.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
ax1.invert_yaxis()
ax1.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax1.set_yticklabels(column_labels)

plt.xlabel("Amino Acid Position")
plt.ylabel("Amino Acid")
plt.title("Amino Acid Fitness")
plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)
plt.tight_layout()
plt.show()


"""
#Plot heatmaps of barcode counts for each amino acid and position
column_labels = [x for x,y in sorted(amino_to_number_dict.items(), key=lambda x: x[1])]
row_labels = sorted(aa_fitness_scores.keys())
ax1 = plt.subplot(211) #this plot will use a linear color scale
heatmap = plt.imshow(aa_fitness_matrix, cmap=plt.cm.YlGnBu, interpolation="nearest")

ax1.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
#ax1.set_yticks(np.arange(barcode_count_aa_heatmap.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
ax1.invert_yaxis()
ax1.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax1.set_yticklabels(column_labels)

plt.xlabel("Amino Acid Position")
plt.ylabel("Amino Acid")
plt.title("Amino Acid Fitness")
plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)
plt.tight_layout()
plt.show()
'''
'''
ax2 = plt.subplot(212) #this plot will use a log color scale
heatmap = plt.imshow(barcode_count_aa_heatmap, norm=LogNorm(vmin=1, vmax=131), cmap=plt.cm.YlGnBu, interpolation="nearest")

ax2.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
ax2.set_yticks(np.arange(barcode_count_aa_heatmap.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))
ax2.invert_yaxis()
ax2.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax2.set_yticklabels(column_labels)

plt.xlabel("Amino Acid Position")
plt.ylabel("Amino Acid")
plt.title("Barcode counts per amino acid at each sequence position (log scale)")
plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)

plt.tight_layout()
plt.show()


#Plot heatmap of deviation from expected proportion of codon at a given aa position
column_labels = [x for x,y in sorted(codon_to_number_dict.items(), key=lambda x: x[1])]
row_labels = sorted_seq_pos

ax = plt.subplot()
heatmap = plt.imshow(barcode_deviation_from_expected, norm=Normalize(vmin=-0.1, vmax=0.1), cmap=plt.cm.seismic, interpolation="nearest")
ax.set_xticks(np.arange(row_labels[0], row_labels[-1], 5))
ax.set_yticks(np.arange(barcode_deviation_from_expected.shape[0]), minor=False)
plt.gca().set_xlim((-0.5, len(row_labels) - 0.5))
plt.gca().set_ylim((-0.5, len(column_labels) - 0.5))

ax.invert_yaxis()
ax.set_xticklabels(scipy.arange(row_labels[0] + 2, row_labels[-1] + 2, 5))
ax.set_yticklabels(column_labels)

plt.tick_params(axis='both', which='major', labelsize=10)

plt.xlabel("Codon Position")
plt.ylabel("Codon")
plt.title("Deviation of proportion barcode at codon-position from expectation with no bias")

plt.colorbar(orientation='vertical', shrink=0.75, pad=0.01)

plt.tight_layout()
plt.show()
'''

'''
#Plot GC count histograms
ax = plt.subplot()
bins = np.arange(0.0,1.1,0.1)
ax.hist(codon_gc_perc, bins=[0, 0.25, 0.75, 1], color='black', label='Codons')
ax.hist(barcode_gc_perc, bins=bins, color='grey', label='Barcodes')
bin_centers = 0.5*(bins[1:]+bins[:-1])
ax.set_xticks(bin_centers)
ax.set_xticklabels(bins)
plt.gca().set_xlim(-0.1,1.1)
plt.title("GC content of barcodes and codons")
plt.xlabel("GC percent")
plt.ylabel("Count")
plt.legend(loc=2)
plt.show()
'''
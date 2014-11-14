import cPickle as pic 
import numpy as np
from create_mat_barcode_fit_dic import *
import sys
import matplotlib.pyplot as plt

# input: dictionary of fitness scores corresponding to a codon/position
# barcode value:[fitness score, standard error, p-value]

# Needs the following pickles in the working folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Needs create_mat_barcode_fit_dict in the working folder.

if __name__ == "__main__":
  filename = sys.argv[1]

  A2N = pic.load(open("aa_to_num_sorted.pkl", "rb"))
  translate = pic.load(open("translate.pkl","rb"))
  codon_dic = pic.load(open('codon_ypos.pkl', "rb"))

  #get matrix 
  fitness_mat = create_mat_barcode_fit_dic(filename)
  iter=0
  big_bc_list = []

  for position in range(len(fitness_mat)):
    wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTATG'

    codon = wtseq[position*3:position*3+3].replace('T', 'U')
    wt_at_pos = translate[codon]
    aa_pos = A2N[wt_at_pos]
    wt_barcodes = fitness_mat[position][aa_pos]
    if len(wt_barcodes) >0:
      for i in wt_barcodes.keys():
        if not np.isnan(wt_barcodes[i]):
          big_bc_list.append(wt_barcodes[i])
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  ax.set_ylabel('Counts')
  ax.set_xlabel('Fitness Score')
  plt.hist(big_bc_list, bins=500)
  plt.savefig('wt_fitness_hist.png')  
  np.savetxt('wt_fitness.txt', big_bc_list)









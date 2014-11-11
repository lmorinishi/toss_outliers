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

# Run
# cd to file + dependencies
# python toss_outliers_v2.py your_fitness_pkl.pkl

# Find the med and MAD of the stop and wt codon fitness values
def find_var_stop_wt(pos_codon_file):
  # Establish fitness matrix, # aa, # codons, wt sequence
  fit_matrix = pos_codon_file
  p, codon = np.shape(fit_matrix)
  wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTATG'

  # Get pkls
  translate = pic.load(open("translate.pkl","rb"))
  A2N = pic.load(open("aa_to_num_sorted.pkl", "rb"))

  # Establish output lists
  stop_med_var = []
  wt_med_var = []

  # Iterate through position indices, find STOP codons at that position, make list of 
  # fitness vals, append to output list (median)
  for pos in range(p):
    dict_barcodes = fit_matrix[pos][A2N['STOP']]
    stop_scores = dict_barcodes.values()
    stop_med_var.append((np.median(stop_scores), get_mad(dict_barcodes)))

    # At position, convert DNA to RNA, and translate to aa
    cod = wtseq[pos*3:pos*3+3].replace('T', 'U')
    wt_at_pos = translate[cod]
    aa_pos = A2N[wt_at_pos]

    # Find WT codons at position, make list of tuples as above
    wt_barcodes = fit_matrix[pos][aa_pos]
    wt_scores = wt_barcodes.values()
    wt_med_var.append((np.median(wt_scores), get_mad(wt_barcodes)))

  # wt_med_list = [item[0] for item in wt_med_var]
  # wt_var_list = [item[1] for item in wt_med_var]
  # stop_med_list = [item[0] for item in stop_med_var]
  # stop_var_list = [item[1] for item in stop_med_var]

  return stop_med_var, wt_med_var


def remove_all_outliers(pos_codon_file, cutoff = 3):
  wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTATG'

  # Load up necessary pkls
  translate = pic.load(open("translate.pkl","rb"))
  codon_dic = pic.load(open('codon_ypos.pkl', "rb"))
  A2N = pic.load(open("aa_to_num_sorted.pkl", "rb"))

  # Establish output lists
  clean_barcodes = []
  dirty_barcodes = []
  
  # Establish fitness matrix, # aa, # codons, wt sequence
  fit_matrix = pos_codon_file
  p, num_codons = np.shape(fit_matrix)
  
  # Iterate through position indices
  for pos in range(p):

    # Convert WT DNA to RNA to AA to #
    cod = wtseq[pos*3:pos*3+3].replace('T', 'U')
    wt_at_pos = translate[cod]
    aa_pos = A2N[wt_at_pos]

    # Find WT barcodes, get MAD
    wt_barcodes = fit_matrix[pos][aa_pos]
    wt_mad = get_mad(wt_barcodes)

    # Iterate through possible codons
    for c in range(num_codons):
      dict_barcodes = fit_matrix[pos][c]

      # Toss outliers if number of barcodes per codon is high enough
      if len(dict_barcodes) > cutoff:
        clean, dirty = get_outliers(dict_barcodes, wt_mad)
        clean_barcodes = clean_barcodes + clean.keys()
        dirty_barcodes = dirty_barcodes + dirty.keys()
      else:
        #implement a good way to deal with low numbers of fitness scores per codon/position
        clean_barcodes = clean_barcodes + dict_barcodes.keys()

  return clean_barcodes, dirty_barcodes


# COMPUTES A HEATMAP WITH THE DIFFERENCE IN VARIANCE  
def var_heatmap(pos_codon_file, cutoff=3):
  wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTATG'

  # Load relevant pkls
  A2N = pic.load(open("aa_to_num_sorted.pkl", "rb"))
  translate = pic.load(open("translate.pkl","rb"))
  codon_dic = pic.load(open('codon_ypos.pkl', "rb"))
  r_codon_dic = dict(zip(codon_dic.values(), codon_dic.keys()))

  # rename position file and get #positions, #codons
  fit_matrix = pos_codon_file
  p, codon = np.shape(fit_matrix)

  # Makes an empty matrix qxq
  q = 21
  fit_aa_before = [[[] for a in range(q)] for pos in range(p)] 
  fit_aa_after = [[[] for a in range(q)] for pos in range(p)]

  #POPULATE MATRIX WITH FITNESS VALUES BY POS AND AA
  for pos in range(p):

    cod = wtseq[pos*3:pos*3+3].replace('T', 'U')
    wt_at_pos = translate[cod]
    aa_pos = A2N[wt_at_pos]
    wt_barcodes = fit_matrix[pos][aa_pos]
    wt_mad = get_mad(wt_barcodes)

    for c in range(codon):
      dict_barcodes = fit_matrix[pos][c]
      if len(dict_barcodes) > cutoff:
        clean, dirty = get_outliers(dict_barcodes, wt_mad)
      else:
        clean = dict_barcodes
        dirty = clean

      am_num = A2N[translate[r_codon_dic[c].replace('T', 'U')]]
      fit_aa_after[pos][am_num] = fit_aa_after[pos][am_num] + clean.values()
      fit_aa_before[pos][am_num] = fit_aa_before[pos][am_num] + dict_barcodes.values()

  #COMPUTE VARIANCES BY POS AND AA
  fit_var_before = np.zeros((p,q))
  fit_var_after = np.zeros((p,q))
  for pos in range(p):
    for a in range(q):
       if len(fit_aa_after[pos][a])>1:
         fit_var_before[pos,a]=np.std(np.array(fit_aa_before[pos][a]))
         fit_var_after[pos,a]=np.std(np.array(fit_aa_after[pos][a]))


  #PLOT HEATMAP 
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1, aspect=(.75*p/q))
  ax.set_yticks(np.array(A2N.values())+.5)
  ax.set_yticklabels(A2N.keys())
  ax.set_xlim([0,p])
  plt.pcolor((fit_var_before.T-fit_var_after.T), cmap='Spectral')
  plt.pcolor((fit_var_after.T), cmap='Spectral')
  plt.colorbar()
  plt.savefig("fitness_variance_heatmap.png")

  return fit_var_before, fit_var_after

# Get MAD of fitness scores. Expects dict {barcode:fitness}
def get_mad(fitness_scores):
  scores = fitness_scores.values()
  scores = np.array(scores)
  #calculate MAD
  med = np.median(scores)
  mad = 1.4826*np.median(abs(np.subtract(scores, med)))
  return mad


def get_outliers(fitness_scores, mad, thresh=2.5):
  
  fit_tup = [(k, fitness_scores[k]) for k in fitness_scores]
  scores = [x[1] for x in fit_tup]

  scores = np.array(scores)
  #calculate MAD
  med = np.median(scores)

  # mad = 1.4826*np.median(abs(np.subtract(scores, med)))

  numerator = abs(scores)- float(med)

  mad_scores = numerator/mad

  final_dict = {}
  outliers = {}

  for i in range(len(mad_scores)):
    key = fit_tup[i][0]
    value = fit_tup[i][1]
  
    if mad_scores[i] < thresh:
      final_dict[key] = value
    else:
      outliers[key] = value

  return final_dict, outliers



if __name__ == "__main__":
  filename = sys.argv[1]

  #open pickle file from other group
  #barcodes = pic.load(open(filename))

  #get matrix 
  fitness_mat = create_mat_barcode_fit_dic(filename)

  #compute heatmap by aa
  before, after = var_heatmap(fitness_mat)
  #remove outliers
  list_barcodes = remove_all_outliers(fitness_mat)

  stop_med_var, wt_med_var = find_var_stop_wt(fitness_mat)

  #^ this is the list we need to give back to the fitness group








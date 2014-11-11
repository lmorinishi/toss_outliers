import cPickle as pic 
import numpy as np
from create_mat_barcode_fit_dic import *
import sys
import matplotlib.pyplot as plt


#input: dictionary of fitness scores corresponding to a codon/position
#key: barcode value: fitness score
#THIS IS FOR A LARGE NUMBER OF BARCODES PER CODON


#sort input barcode dictionary to find barcodes corresponding to a codon/position

def remove_all_outliers(pos_codon_file, cutoff = 6):

	clean_barcodes = []
	dirty_barcodes = []

	fit_matrix = pos_codon_file
	p, codon = np.shape(fit_matrix)
	#check position indices
	for pos in range(p):
		for c in range(codon):
			dict_barcodes = fit_matrix[pos][c]
			if len(dict_barcodes) > cutoff:
				clean, dirty = get_outliers(dict_barcodes)
				clean_barcodes = clean_barcodes + clean.keys()
				dirty_barcodes = dirty_barcodes + dirty.keys()
			else:
				#implement a good way to deal with low numbers of fitness scores per codon/position
				clean_barcodes = clean_barcodes + dict_barcodes.keys()
	return clean_barcodes, dirty_barcodes


#COMPUTES A HEATMAP WITH THE DIFFERENCE IN VARIANCE BETWEEN 
def var_heatmap(pos_codon_file, cutoff=6):

  A2N = pic.load(open("aa_to_num_sorted.pkl", "rb"))
  translate = pic.load(open("translate.pkl","rb"))
  codon_dic = pic.load(open('codon_ypos.pkl', "rb"))
  r_codon_dic = dict(zip(codon_dic.values(), codon_dic.keys()))

  fit_matrix = pos_codon_file
  p, codon = np.shape(fit_matrix)
  q = 21
  print p,q
  fit_aa_before = [[[] for a in range(q)] for pos in range(p)] 
  fit_aa_after = [[[] for a in range(q)] for pos in range(p)]

  #POPULATE MATRIX WITH FITNESS VALUES BY POS AND AA
  for pos in range(p):
     for c in range(codon):
        dict_barcodes = fit_matrix[pos][c]
        if len(dict_barcodes) > cutoff:
           clean, dirty = get_outliers(dict_barcodes)
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
  plt.pcolormesh((fit_var_before.T-fit_var_after.T), cmap='Spectral')
  plt.colorbar()
  plt.savefig("fitness_variance_heatmap.png")

  return fit_var_before, fit_var_after

def get_outliers(fitness_scores, thresh=2.5):
      
      fit_tup = [(k, fitness_scores[k]) for k in fitness_scores]
      scores = [x[1] for x in fit_tup]

      scores = np.array(scores)
      #calculate MAD
      med = np.median(scores)


      mad = 1.4826*np.median(abs(np.subtract(scores, med)))

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

	#^ this is the list we need to give back to the fitness group








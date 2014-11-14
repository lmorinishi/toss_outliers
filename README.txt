# input: dictionary of fitness scores corresponding to a codon/position
# barcode value:[fitness score, standard error, p-value]

# Needs the following pickles in the working folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Input: dictionary of fitness scores corresponding to a codon or mutation
	# Either:
	# {barcode:[fitness, stde, pval, rval]}
	# {(pos,codon):[fitness, stde, pval, rval]}

# Run
# cd to directory
# python toss_outliers_v3.py pickle.pkl sort_type cutoff
# where sort_type should be either 'barcode' or 'codon'
# where cutoff should be an int, probably 4 for barcode and 0 for codon


# Test fitness of wt barcodes
# python wt_hist.py your_fitness_pkl.pkl
# Saves a histogram wt_fitness_hist.png and list of fitness values wt_fitness.txt
# python distribution-check.py wt_fitness.txt
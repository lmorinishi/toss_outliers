# Needs the following pickles in the pkl folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Needs pkls to be analyzed in fitness folder.

# Input: dictionary of fitness scores corresponding to a codon or mutation
	# Either:
	# {barcode:[fitness, stde, pval, rval]}
	# {(pos,codon):[fitness, stde, pval, rval]}

# Run
# cd to directory
# ./run.sh
# Input barcode cutoff (recommend 4) and codon cutoff (rec 1)
# Outputs /output/blah_clean_list.pkl
	# {clean_barcodes:[barcodes], dirty_barcodes[barcodes]} OR
	# {clean_barcodes:[(pos,aa)], dirty_barcodes[(pos,aa)]}
# Also writes number of “clean” and “dirty (outlier)” barcodes in list_counts.txt


# Test fitness of wt barcodes
# python wt_hist.py your_fitness_pkl.pkl
# Saves a histogram wt_fitness_hist.png and list of fitness values wt_fitness.txt
# python distribution-check.py wt_fitness.txt
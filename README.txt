# input: dictionary of fitness scores corresponding to a codon/position
# barcode value:[fitness score, standard error, p-value]

# Needs the following pickles in the working folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Needs create_mat_barcode_fit_dict in the working folder.

# Run
# cd to file + dependencies
# python toss_outliers_v1.py your_fitness_pkl.pkl


# Test fitness of wt barcodes
# python wt_hist.py your_fitness_pkl.pkl
# Saves a histogram wt_fitness_hist.png and list of fitness values wt_fitness.txt
# python distribution-check.py wt_fitness.txt
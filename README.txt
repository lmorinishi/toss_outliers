# input: dictionary of fitness scores corresponding to a codon/position
# barcode value:[fitness score, standard error, p-value]

# Needs the following pickles in the working folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Run
# cd to file + dependencies
# python toss_outliers_v2.py your_fitness_pkl.pkl
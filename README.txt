# Needs the following pickles in the pkl folder.
# translate.pkl, codon_ypos.pkl, aa_to_num_sorted.pkl, allele_dic_with_WT.pkl

# Needs pkls to be analyzed in fitness folder.

# Input: dictionary of fitness scores in the following format
	# {barcode:[fitness, stde, pval, rval]}

# Run
# cd to directory
# run python toss_outliers_aa.py <input_pkl> <cutoff> <output_directory_with_slash>
#EXAMPLE: python toss_outliers_aa.py fitness/Barcode_FitScore_Caffiene_day1.pkl 5 /Users/ipqb/Documents/Fall-2014/PUBS/barcode-reproducibility/

# Outputs /output/blah_clean_list.pkl
	# {clean_barcodes:[barcodes], dirty_barcodes:[barcodes]} 
#Also outputs some plots to/heatmap/
# Also writes number of “clean” and “dirty (outlier)” barcodes in list_counts.txt

# Test fitness of wt barcodes
# python wt_hist.py your_fitness_pkl.pkl
# Saves a histogram wt_fitness_hist.png and list of fitness values wt_fitness.txt
# python distribution-check.py wt_fitness.txt
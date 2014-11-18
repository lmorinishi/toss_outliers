import cPickle as pic 
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Input: dictionary of fitness scores corresponding to a codon or mutation
	# Either:
	# {barcode:[fitness, stde, pval, rval]}
	# {(pos,codon):[fitness, stde, pval, rval]}

# Run
# cd to directory
# python toss_outliers_aa.py pickle.pkl sort_type cutoff
# where sort_type should be either 'barcode' or 'codon'
# where cutoff should be an int, probably 4 for barcode and 0 for codon

# CREATES A MATRIX OF DICTIONARIES, VALUE=FITNESS, SIZE AND KEYS VARY BY SORT_TYPE
A2N = pic.load(open('pkl/aa_to_num_sorted.pkl',"rb"))
bc_mut_dict = pic.load(open('pkl/allele_dic_with_WT.pkl',"rb"))
translate = pic.load(open('pkl/translate.pkl',"rb"))
codon_ypos = pic.load(open('pkl/codon_ypos.pkl', 'rb'))
wtseq='ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTTGA'

#two outputs, codon matrix and aa matrix
def create_fitness_mat(filename, sort_type='barcode'):

	# Calling all pkls
	fitness_dict = pic.load(open(filename, 'rb'))

	# Establish number of aa positions
	p = 77
	q = 64
	a=21

	wt_barcodes=[]

	# Make Npos x Ncodon matrix of dictionaries, {barcode:fitness}

	fitness_mat = [[{} for i in range(q)] for j in range(p)]
	fitness_aa= [[{} for i in range(a)] for j in range(p)]

	for barcode in fitness_dict.keys():
		pos, codon = bc_mut_dict[barcode]
		if codon == 'WT':
			wt_barcodes.append(barcode)

		else:
			rna = codon.replace('T','U')
			aa = translate[rna]
			if not np.isnan(fitness_dict[barcode][0]):
				fitness_mat[pos-1][codon_ypos[codon]][barcode] = fitness_dict[barcode][0]
				fitness_aa[pos-1][A2N[aa]][barcode]= fitness_dict[barcode][0]

				
			else:
				print fitness_dict[barcode] , 

	
	#go back and add in the WT fitness scores into the appropriate positions in the matrix

	p=1
	for i in range(0, len(wtseq), 3):
		cod=wtseq[i:i+3]
		rna=cod.replace('T', 'U')
		aa=translate[rna]
		for b in wt_barcodes:
			fitness_mat[p-1][codon_ypos[cod]][b]=fitness_dict[b][0]
			fitness_aa[p-1][A2N[aa]][b]=fitness_dict[b][0]

		p+=1
	return fitness_mat, fitness_aa


def remove_all_outliers(pos_aa_file, cutoff = 5):
	clean_barcodes = []
	dirty_barcodes = []

	fit_matrix = pos_aa_file
	p, 	q = np.shape(fit_matrix)

	before_mat = [[[] for i in range(q)] for j in range(p)]
	after_mat = [[[] for i in range(q)] for j in range(p)]


	#print p, aa
	#^ check position indices

	for pos in range(p):
		for a in range(q):
			dict_barcodes = fit_matrix[pos][a]
			if len(dict_barcodes) > cutoff:
				clean, dirty = get_outliers(dict_barcodes)
				clean_barcodes = clean_barcodes + clean.keys()
				dirty_barcodes = dirty_barcodes + dirty.keys()
			else:
				#no filtering is done, all barcodes are kept
				clean=dict_barcodes
				dirty={}
				clean_barcodes = clean_barcodes + dict_barcodes.keys()

			for fitness in dict_barcodes.values():
				before_mat[pos][a].append(fitness)

			for barcode in clean.keys():
				after_mat[pos][a].append(clean[barcode])


	return clean_barcodes, dirty_barcodes, before_mat, after_mat


#use an MAD threshold of 2.. filters out a lot of barcodes
#recommended MAD threshold is between 2.0 and 3.0. 

def plot_results(before_mat, after_mat, output_dir, output_file):

	#get variance before and variance after
	p, q = np.shape(before_mat)

	fit_var_before = np.zeros((p,q))
	fit_var_after = np.zeros((p,q))

	num_before=np.add(np.zeros((p,q)),1)
	num_after=np.add(np.zeros((p,q)),1)

	for pos in range(p):
		for a in range(q):
			#variance and counts
			if len(before_mat[pos][a]) >= 1:
				fit_var_before[pos][a]=np.std(np.array(before_mat[pos][a]))
				num_before[pos][a]=len(before_mat[pos][a])

			if len(after_mat[pos][a]) >= 1:
				fit_var_after[pos][a]=np.std(np.array(after_mat[pos][a]))
				num_after[pos][a]=len(after_mat[pos][a])




	minv = min([np.amin(fit_var_before)] + [np.amin(fit_var_after)])
	maxv = max([np.amax(fit_var_before)] + [np.amax(fit_var_after)])
	# w = [x/10.0 for x in range(-1,4)]

	# IMPORTANT FOR ACCURATE TICKS ON HEATMAP
	n = [.5]+list(np.arange(4.5,76.5,5))
	m = [1]+list(np.arange(5,77,5))

	fig = plt.figure()
	fig.subplots_adjust(hspace=.3)

	# PLOT HEATMAPS
	p = 77
	q = 21

	# Plot heatmaps
	ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
	ax.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax.tick_params(axis='both', which='minor', labelsize=8)
	ax.tick_params(axis='both', which='major', labelsize=5)
	ax.set_xticks(n, minor=True)
	ax.set_xticklabels(m, minor=True)
	ax.set_yticks(np.array(A2N.values())+.5)
	ax.set_yticklabels(A2N.keys())
	ax.set_xlim([0,p])
	ax.set_title('Before tossing outliers', fontsize=10)
	ax.set_aspect('equal')
	plt.pcolormesh((fit_var_before.T), cmap='Spectral', vmin=minv, vmax=maxv)
	cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
	cbar.ax.tick_params(labelsize=10) 

	ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
	ax2.tick_params(axis='both', which='minor', labelsize=8)
	ax2.tick_params(axis='both', which='major', labelsize=5)
	ax2.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax2.set_xticks(n, minor=True)
	ax2.set_xticklabels(m, minor=True)
	ax2.set_yticks(np.array(A2N.values())+.5)
	ax2.set_yticklabels(A2N.keys())
	ax2.set_xlim([0,p])
	ax2.set_title('After tossing outliers', fontsize=10)
	ax2.set_aspect('equal')
	plt.pcolormesh((fit_var_after.T), cmap='Spectral', vmin=minv, vmax=maxv)
	cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
	cbar2.ax.tick_params(labelsize=10) 

	ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
	ax3.tick_params(axis='both', which='minor', labelsize=8)
	ax3.tick_params(axis='both', which='major', labelsize=5)
	ax3.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax3.set_xticks(n, minor=True)
	ax3.set_xticklabels(m, minor=True)
	ax3.set_yticks(np.array(A2N.values())+.5)
	ax3.set_yticklabels(A2N.keys())
	ax3.set_xlim([0,p])
	ax3.set_title('Subtraction of before-after', fontsize=10)
	ax3.set_aspect('equal')
	plt.pcolormesh((fit_var_before.T-fit_var_after.T), cmap='Spectral')
	cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
	cbar3.ax.tick_params(labelsize=10) 
	
	heatmaps_d = output_dir+'/heatmaps/'
	if not os.path.exists(heatmaps_d):
		os.makedirs(heatmaps_d)

	plt.savefig(heatmaps_d+output_file+'_fitness_variance_by_aa.png')
	plt.close()




	fig=plt.figure()
	fig.subplots_adjust(hspace=.3)
	#plot number of barcodes before and after tossing outliers
	ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
	ax.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax.tick_params(axis='both', which='minor', labelsize=8)
	ax.tick_params(axis='both', which='major', labelsize=5)
	ax.set_xticks(n, minor=True)
	ax.set_xticklabels(m, minor=True)
	ax.set_yticks(np.array(A2N.values())+.5)
	ax.set_yticklabels(A2N.keys())
	ax.set_xlim([0,p])
	ax.set_title('Before tossing outliers', fontsize=10)
	ax.set_aspect('equal')
	plt.pcolormesh((np.log2(num_before).T), cmap='Spectral')
	cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
	cbar.ax.tick_params(labelsize=10)

	ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
	ax2.tick_params(axis='both', which='minor', labelsize=8)
	ax2.tick_params(axis='both', which='major', labelsize=5)
	ax2.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax2.set_xticks(n, minor=True)
	ax2.set_xticklabels(m, minor=True)
	ax2.set_yticks(np.array(A2N.values())+.5)
	ax2.set_yticklabels(A2N.keys())
	ax2.set_xlim([0,p])
	ax2.set_title('After tossing outliers', fontsize=10)
	ax2.set_aspect('equal')
	plt.pcolormesh((np.log2(num_after).T), cmap='Spectral')
	cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
	cbar2.ax.tick_params(labelsize=10)


	ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
	ax3.tick_params(axis='both', which='minor', labelsize=8)
	ax3.tick_params(axis='both', which='major', labelsize=5)
	ax3.set_xticklabels('')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	ax3.set_xticks(n, minor=True)
	ax3.set_xticklabels(m, minor=True)
	ax3.set_yticks(np.array(A2N.values())+.5)
	ax3.set_yticklabels(A2N.keys())
	ax3.set_xlim([0,p])
	ax3.set_title('Subtraction of before-after', fontsize=10)
	ax3.set_aspect('equal')

	mat=num_before.T-num_after.T
	mat[mat==0] = 1
	plt.pcolormesh((np.log2(mat)), cmap='Spectral')
	cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
	cbar3.ax.tick_params(labelsize=10) 

	plt.savefig(heatmaps_d+output_file+'_counts_by_aa.png')

	plt.close()


def get_outliers(fitness_dict, thresh=2.0):
			
	fit_tup = [(k, fitness_dict[k]) for k in fitness_dict]
	scores = np.array([x[1] for x in fit_tup])

	#calculate MAD
	med = np.median(scores)
	mad = 1.4826*np.median(abs(np.subtract(scores, med)))
	numerator = abs(scores- float(med))
	mad_scores = numerator/mad

	final_dict = {}
	outliers = {}

	for i in range(len(mad_scores)):
		key = fit_tup[i][0]
		fitness = fit_tup[i][1]
			
		if mad_scores[i] < thresh:
			final_dict[key] = fitness
		else:
			outliers[key] = fitness
	return final_dict, outliers

# COMPUTES A HEATMAP OF THE DIFFERENCE IN VARIANCE 
# def var_heatmap(fitness_mat,  output_dir, cutoff=5):

# 	A2N = pic.load(open("pkl/aa_to_num_sorted.pkl", "rb"))
# 	codon_ypos = pic.load(open("pkl/codon_ypos.pkl", "rb"))
# 	translate = pic.load(open("pkl/translate.pkl","rb"))
# 	codon_dic = pic.load(open('pkl/codon_ypos.pkl', "rb"))
# 	pathname = os.path.split(filename)

# 	# Determine shape of matrix & make empty lists
# 	p, q = np.shape(fitness_mat)
# 	print p, q

# 	before_mat = [[[] for i in range(q)] for j in range(p)]
# 	after_mat = [[[] for i in range(q)] for j in range(p)]
# 	clean_list = []
# 	dirty_list = []



# 	# Iterate through matrix, get outliers for each
# 	# separate into clean and dirty dictionaries for heatmap
# 	# append to lists for ultimate pkl output
# 	for pos in range(p):
# 		for c in range(q):
# 			dict_barcodes = fitness_mat[pos][c]
# 			if len(dict_barcodes) >= cutoff:
# 				clean, dirty = get_outliers(dict_barcodes)
# 			else:
# 				clean = dict_barcodes
# 				dirty = {}

# 			for barcode in clean.keys():
# 				clean_list.append(barcode)
# 				after_mat[pos][c].append(clean[barcode])
# 			for fitness in dict_barcodes.values():
# 				before_mat[pos][c].append(fitness)
# 			for barcode in dirty.keys():
# 				dirty_list.append(barcode)

# 	#COMPUTE VARIANCES BY POS AND AA/CODON
# 	fit_var_before = np.zeros((p,q))
# 	fit_var_after = np.zeros((p,q))
# 	num_before=np.zeros((p,q))
# 	num_after=np.zeros((p,q))

# 	for pos in range(p):
# 		for c in range(q):
# 			if len(before_mat[pos][c]) != 0:
# 				num_before[pos][c] = np.log2(len(before_mat[pos][c]))
# 			else:
# 				num_before[pos][c] = 0
# 			if len(after_mat[pos][c])!= 0:
# 				num_after[pos][c] = np.log2(len(after_mat[pos][c]))
# 			else:
# 				num_after[pos][c] = 0

# 			if len(after_mat[pos][c]) > 1:
# 				fit_var_after[pos,c]=np.std(np.array(after_mat[pos][c]))
# 			if len(before_mat[pos][c]) > 1:
# 				fit_var_before[pos,c]=np.std(np.array(before_mat[pos][c]))	

# 	# ESTABLISH MIN/MAX TO NORMALIZE HEATMAPS
# 	minv = min([np.amin(fit_var_before)] + [np.amin(fit_var_after)])
# 	maxv = max([np.amax(fit_var_before)] + [np.amax(fit_var_after)])
# 	# w = [x/10.0 for x in range(-1,4)]

# 	# IMPORTANT FOR ACCURATE TICKS ON HEATMAP
# 	n = [.5]+list(np.arange(4.5,76.5,5))
# 	m = [1]+list(np.arange(5,77,5))

# 	fig = plt.figure()
# 	fig.subplots_adjust(hspace=.3)

# 	# PLOT HEATMAPS
# 	p = 77
# 	q = 21

# 	# Plot heatmaps
# 	ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
# 	ax.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax.tick_params(axis='both', which='minor', labelsize=8)
# 	ax.tick_params(axis='both', which='major', labelsize=5)
# 	ax.set_xticks(n, minor=True)
# 	ax.set_xticklabels(m, minor=True)
# 	ax.set_yticks(np.array(A2N.values())+.5)
# 	ax.set_yticklabels(A2N.keys())
# 	ax.set_xlim([0,p])
# 	ax.set_title('Before tossing outliers', fontsize=10)
# 	ax.set_aspect('equal')
# 	plt.pcolormesh((fit_var_before.T), cmap='Spectral', vmin=minv, vmax=maxv)
# 	cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
# 	cbar.ax.tick_params(labelsize=10) 

# 	ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
# 	ax2.tick_params(axis='both', which='minor', labelsize=8)
# 	ax2.tick_params(axis='both', which='major', labelsize=5)
# 	ax2.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax2.set_xticks(n, minor=True)
# 	ax2.set_xticklabels(m, minor=True)
# 	ax2.set_yticks(np.array(A2N.values())+.5)
# 	ax2.set_yticklabels(A2N.keys())
# 	ax2.set_xlim([0,p])
# 	ax2.set_title('After tossing outliers', fontsize=10)
# 	ax2.set_aspect('equal')
# 	plt.pcolormesh((fit_var_after.T), cmap='Spectral', vmin=minv, vmax=maxv)
# 	cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
# 	cbar2.ax.tick_params(labelsize=10) 

# 	ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
# 	ax3.tick_params(axis='both', which='minor', labelsize=8)
# 	ax3.tick_params(axis='both', which='major', labelsize=5)
# 	ax3.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax3.set_xticks(n, minor=True)
# 	ax3.set_xticklabels(m, minor=True)
# 	ax3.set_yticks(np.array(A2N.values())+.5)
# 	ax3.set_yticklabels(A2N.keys())
# 	ax3.set_xlim([0,p])
# 	ax3.set_title('Subtraction of before-after', fontsize=10)
# 	ax3.set_aspect('equal')
# 	plt.pcolormesh((fit_var_before.T-fit_var_after.T), cmap='Spectral')
# 	cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
# 	cbar3.ax.tick_params(labelsize=10) 
	
# 	heatmaps_d = output_dir+'/heatmaps/'
# 	if not os.path.exists(heatmaps_d):
# 		os.makedirs(heatmaps_d)

# 	plt.savefig(heatmaps_d+pathname[1]+'_fitness_variance_by_aa.png')
# 	plt.close()


# 	fig=plt.figure()
# 	fig.subplots_adjust(hspace=.3)
# 	#plot number of barcodes before and after tossing outliers
# 	ax = fig.add_subplot(3,1,1, aspect=(.75*p/q))
# 	ax.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax.tick_params(axis='both', which='minor', labelsize=8)
# 	ax.tick_params(axis='both', which='major', labelsize=5)
# 	ax.set_xticks(n, minor=True)
# 	ax.set_xticklabels(m, minor=True)
# 	ax.set_yticks(np.array(A2N.values())+.5)
# 	ax.set_yticklabels(A2N.keys())
# 	ax.set_xlim([0,p])
# 	ax.set_title('Before tossing outliers', fontsize=10)
# 	ax.set_aspect('equal')
# 	plt.pcolormesh((num_before.T), cmap='Spectral')
# 	cbar = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
# 	cbar.ax.tick_params(labelsize=10)

# 	ax2 = fig.add_subplot(3,1,2, aspect=(.75*p/q))
# 	ax2.tick_params(axis='both', which='minor', labelsize=8)
# 	ax2.tick_params(axis='both', which='major', labelsize=5)
# 	ax2.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax2.set_xticks(n, minor=True)
# 	ax2.set_xticklabels(m, minor=True)
# 	ax2.set_yticks(np.array(A2N.values())+.5)
# 	ax2.set_yticklabels(A2N.keys())
# 	ax2.set_xlim([0,p])
# 	ax2.set_title('After tossing outliers', fontsize=10)
# 	ax2.set_aspect('equal')
# 	plt.pcolormesh((num_after.T), cmap='Spectral')
# 	cbar2 = plt.colorbar(orientation='vertical', pad=0.01)# ticks=v)
# 	cbar2.ax.tick_params(labelsize=10)


# 	ax3 = fig.add_subplot(3,1,3, aspect=(.75*p/q))
# 	ax3.tick_params(axis='both', which='minor', labelsize=8)
# 	ax3.tick_params(axis='both', which='major', labelsize=5)
# 	ax3.set_xticklabels('')
# 	plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 	ax3.set_xticks(n, minor=True)
# 	ax3.set_xticklabels(m, minor=True)
# 	ax3.set_yticks(np.array(A2N.values())+.5)
# 	ax3.set_yticklabels(A2N.keys())
# 	ax3.set_xlim([0,p])
# 	ax3.set_title('Subtraction of before-after', fontsize=10)
# 	ax3.set_aspect('equal')
# 	plt.pcolormesh((num_before.T-num_after.T), cmap='Spectral')
# 	cbar3 = plt.colorbar(orientation='vertical', pad=0.01)#, ticks=w)
# 	cbar3.ax.tick_params(labelsize=10) 

# 	plt.savefig(heatmaps_d+pathname[1]+'_counts_by_aa.png')


# 	plt.close()

# 	return before_mat, after_mat, clean_list, dirty_list


if __name__ == "__main__":
	filename = sys.argv[1]
	#sort_type = sys.argv[2]
	cutoff = int(sys.argv[2])
	output_f=sys.argv[3]

	# make matrix of dictionaries, with fitness as values
	fitness_mat, fitness_aa = create_fitness_mat(filename)

	# compute and output heatmap by aa


	clean_barcodes, dirty_barcodes, before_mat, after_mat = remove_all_outliers(fitness_aa, cutoff)
	plot_results(before_mat, after_mat, output_f, os.path.split(filename)[1])


	# PRINT NUMBER OF CLEAN AND DIRTY BARCODES FOUND, UNCOMMENT BELOW TO ALSO
	# PRINT NUMBER OF FITNESS SCORES THROWN OUT
	# fitness_dict = pic.load(open(filename, 'rb'))
	# x = len(fitness_dict) - len(clean_list) - len(dirty_list)
	print len(clean_barcodes), '\tclean barcodes'
	print len(dirty_barcodes), '\tdirty barcodes'
	# print x, '\tWTs thrown out'

	# OUTPUT DICT
	pathname = os.path.split(filename)

	output_dict = {'clean_barcodes':clean_barcodes, 'dirty_barcodes':dirty_barcodes}

	out_d = output_f+'/output/'
	if not os.path.exists(out_d):
		os.makedirs(out_d)
	pic.dump(output_dict, open(out_d + pathname[1] + '_clean_list.pkl', 'w'), protocol=2)

	with open(out_d+'list_counts.txt', 'a') as tossfile:
		# txt file: pkl_file clean_counts dirty_counts
		newline = [pathname[1], str(len(clean_barcodes)), str(len(dirty_barcodes))]
		catline = '\t'.join(newline)
		tossfile.write(catline+'\n')


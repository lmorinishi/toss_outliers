import cPickle as pic
import numpy as np
import matplotlib as plt
import pylab as pl

#Creates a matrix of position vs codon with a list of tuples containing (barcode, fitness score)
def create_mat_barcode_fit_dic(pkl_to_use_path):
        
    # primer_pkl = "./pickles_files/ATTGGC.pkl"
    #pkl_to_use_path = "./pickles_files/ACATCG.pkl"
    barcode_mat = pic.load(open(pkl_to_use_path,"rb")) 
    aa = pic.load(open("aa_to_num_sorted.pkl", "rb"))
    barcode_mut_dict = pic.load(open("allele_dic_with_WT.pkl","rb"))
    translate = pic.load(open("translate.pkl","rb"))
    codon_dic = pic.load(open('codon_ypos.pkl', "rb"))

    wtseq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTATG'


    # rows = len(codon_dic)   #number of codons                           
    # columns = 77            #length of ubiquitin            
    # mt_counts = np.ones([rows,columns])

    rows_aa = 21
    columns_aa = 77
    barcode_pos_aa_val_dic = {}
    codon_matrix = [ [ {} for i in range(rows_aa) ] for j in range(columns_aa) ]


    #codon dic maps codon -> position
    #key is barcode, val is fitness score
    for barcode, fitness in barcode_mat.iteritems():
        if barcode in barcode_mut_dict:
            #return value from barcode_mut_dic is {position of mutation, codon mutant}
            pos = barcode_mut_dict[barcode][0]

            if not(barcode_mut_dict[barcode][1] == 'WT'):
                cod = barcode_mut_dict[barcode][1]
                aa_for_cod = translate[cod.replace('T', 'U')]
                codon_matrix[pos-1][aa[aa_for_cod]][barcode] = fitness
            else: 
                cod = wtseq[pos*3:pos*3+3].replace('T', 'U')
                aa_for_cod = translate[cod]
                codon_matrix[pos][aa[aa_for_cod]][barcode] = fitness
            
    return codon_matrix





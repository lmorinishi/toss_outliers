#!/bin/bash

echo `pwd`

echo -e 'Barcode cutoff: \c '
read bc_cutoff
echo -e 'Codon cutoff: \c '
read codon_cutoff

for bc_pkl in $( ls fitness/Barcode*day*.pkl ); do
	echo $bc_pkl
	echo `python toss_outliers_v3.py $bc_pkl 'barcode' $bc_cutoff`
done

for codon_pkl in $( ls fitness/Codon*day*.pkl ); do
	echo $codon_pkl
	echo `python toss_outliers_v3.py $codon_pkl 'codon' $codon_cutoff`
done
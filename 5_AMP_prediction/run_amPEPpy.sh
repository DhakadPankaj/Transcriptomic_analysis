#! /bin/bash
sconda braker3
ampep train --test-trees --min-trees 23 --max-trees 175 -p AMP.tr.fa -n DECOY.tr.fa --seed 100 > oob_error_results.tsv

sort -k2,2n oob_error_results.tsv 

ampep train --feature-importance --num-trees 106 -p AMP.tr.fa -n DECOY.tr.fa --seed 100 > feature.importances.dropcol.oob.csv

#nano features.remove.txt
ampep train -p AMP.tr.fa -n DECOY.tr.fa --drop-features features.remove.txt --seed 100

ampep predict -m amPEP.model -i NA_genes.fasta --drop-features features.remove.txt -o results.tsv --seed 100
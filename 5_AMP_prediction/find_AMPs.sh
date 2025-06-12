#! /bin/bash


# This script is used to get fasta sequence of Unknown short sequences (AMPs)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

for sp in "Hcam" "Hconf" "Sdef"
do
tail -n +2 DESeq2_genes/${sp}_Dmel_genes.tsv| awk '$4=="NA" && $9>=1'|cut -f2 > DESeq2_genes/${sp}_NA_genes.tmp

# get fasta sequence of Unknown short sequences (AMPs), length =< 200
seqtk subseq <(seqkit seq -g -M 200 DESeq2_genes/${sp}_UP.fasta) DESeq2_genes/${sp}_NA_genes.tmp > DESeq2_genes/${sp}_NA_genes.fasta

rm DESeq2_genes/${sp}_NA_genes.tmp

done
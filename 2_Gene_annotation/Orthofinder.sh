#! /bin/bash

while [ "$1" != "" ] ; do
    case $1 in
        -isoforms )           shift
                                isoforms=$1
                                ;;
        -orthofinder )         shift
                                orthofinder=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done

if [ "$isoforms" = "true" ]; then
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agat
mkdir -p orthofinder
# get longest isoform per gene
while IFS=$'\t' read species sp genome; do
    #agat_sp_keep_longest_isoform.pl -gff annot/${sp}/braker.gff3 -o annot/${sp}/${sp}_longest_isoforms.gff
    #gffread -y - -C -S -g genomes/${genome} annot/${sp}/${sp}_longest_isoforms.gff | fasta_formatter|sed "/^>/ s/>/>${sp}_/" > orthofinder/${sp}.fa
    for treat in "infected" "naive" ; do
    agat_sp_keep_longest_isoform.pl -gff annot/${sp}_${treat}/braker.gff3 -o annot/${sp}_${treat}/${sp}_${treat}_longest_isoforms.gff

    letter=$(echo $treat | cut -c1|tr [a-z] [A-Z])
    gffread -y - -C -S -g genomes/${genome} annot/${sp}_${treat}/${sp}_${treat}_longest_isoforms.gff | fasta_formatter|sed "/^>/ s/>/>${sp}_${letter}_/" > orthofinder/${sp}_${letter}.fa
    done
done < species.txt

fi


if [ "$orthofinder" = "true" ]; then
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
# Run Orthofinder
orthofinder -f orthofinder/ -t 60 -a 12 -M msa -S blast -T iqtree -A mafft -X

fi
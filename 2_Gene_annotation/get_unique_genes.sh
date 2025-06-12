#! /bin/bash

# Description of files
# A and B : Genes with no genomic overlap between infected and naive AND genes not assigned to any orthogroup 
# C : Genes in both infected and naive but not in any orthogroup
# D : Genes in infected only and in a orthogroup
# E : Genes in naive only and in a orthogroup


#cat Phylogenetic_Hierarchical_Orthogroups/N5.tsv |cut -f1,6,7|awk -F"\t" '$2!="" && $3!=""'|tail -n +2|cut -f1,2|awk -F'\t' '{split($2, genes, ", "); for (i in genes) { print $1 "\t" genes[i] }}'|sed 's/\.t[0-9]*//'


get_genes(){
    SPECIES=$1
    OUTPUT_DIR="annot"

    # Genes with no genomic overlap between infected and naive
    cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_1_0_id_list.txt |cut -d "|" -f1|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_I_/" > $OUTPUT_DIR/${SPECIES}_I.txt
    cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_0_1_id_list.txt |cut -d "|" -f2|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_N_/" > $OUTPUT_DIR/${SPECIES}_N.txt

    # Genes in Unassigned orthogroups
    col_I=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups_UnassignedGenes.tsv|tr "\t" "\n"|grep -w -n "${SPECIES}_I"|cut -d ":" -f1)
    col_N=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups_UnassignedGenes.tsv|tr "\t" "\n"|grep -w -n "${SPECIES}_N"|cut -d ":" -f1)

    tail -n +2 orthofinder/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups_UnassignedGenes.tsv| cut -f${col_I} |grep -v "^$"|sed 's/\.t[0-9]*//' > $OUTPUT_DIR/${SPECIES}_I_Unassigned.txt
    tail -n +2 orthofinder/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups_UnassignedGenes.tsv| cut -f${col_N} |grep -v "^$"|sed 's/\.t[0-9]*//' > $OUTPUT_DIR/${SPECIES}_N_Unassigned.txt

    # A and B : Genes with no genomic overlap between infected and naive AND genes not assigned to any orthogroup
    grep -Fwf $OUTPUT_DIR/${SPECIES}_I.txt $OUTPUT_DIR/${SPECIES}_I_Unassigned.txt > $OUTPUT_DIR/${SPECIES}_I_A.txt
    grep -Fwf $OUTPUT_DIR/${SPECIES}_N.txt $OUTPUT_DIR/${SPECIES}_N_Unassigned.txt > $OUTPUT_DIR/${SPECIES}_N_B.txt
    
    # C : Genes in both infected and naive but not in any orthogroup
    grep -Fwf $OUTPUT_DIR/${SPECIES}_I_Unassigned.txt <(cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_[1-9]_[1-9]_id_list.txt |cut -d "|" -f1|tr "," "\n"|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_I_/") > $OUTPUT_DIR/${SPECIES}_I_C.txt
    grep -Fwf $OUTPUT_DIR/${SPECIES}_N_Unassigned.txt <(cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_[1-9]_[1-9]_id_list.txt |cut -d "|" -f2|tr "," "\n"|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_N_/") > $OUTPUT_DIR/${SPECIES}_N_C.txt
    
    # D : Genes in infected only and in a orthogroup
    grep -Fwvf $OUTPUT_DIR/${SPECIES}_I_Unassigned.txt $OUTPUT_DIR/${SPECIES}_I.txt |grep -Fwvf <(tail -n +2 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|awk -F'\t' -v I=$col_I -v N=$col_N '$N != "" && $I != ""'|cut -f${col_I}|tr "," "\n"|sed 's/\s//g'|grep -v '^$'|sed 's/\.t[0-9]*//') > $OUTPUT_DIR/${SPECIES}_I_D.txt
    # E : Genes in naive only and in a orthogroup
    grep -Fwvf $OUTPUT_DIR/${SPECIES}_N_Unassigned.txt $OUTPUT_DIR/${SPECIES}_N.txt |grep -Fwvf <(tail -n +2 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|awk -F'\t' -v I=$col_I -v N=$col_N '$N != "" && $I != ""'|cut -f${col_N}|tr "," "\n"|sed 's/\s//g'|grep -v '^$'|sed 's/\.t[0-9]*//') > $OUTPUT_DIR/${SPECIES}_N_E.txt

    
    
}

export -f get_genes
mkdir -p annot
parallel -j 3 --colsep "\t" get_genes {2} :::: species.txt

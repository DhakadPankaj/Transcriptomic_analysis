#! /bin/bash

# This script is used to get the genes that are in the Venn diagram
# Compare infected annotations with uninfected annotations


get_genes(){
    SPECIES=$1
    OUTPUT_DIR="annot"
    hog_dir="orthofinder/OrthoFinder/Results_Apr05/Phylogenetic_Hierarchical_Orthogroups"

    # Genes with no genomic overlap between infected and naive
    cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_1_0_id_list.txt |cut -d "|" -f1|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_I_/" > $OUTPUT_DIR/${SPECIES}_I.txt
    cat $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_0_1_id_list.txt |cut -d "|" -f2|sed 's/\s//g'|grep -v '^$'| sed "s/^/${SPECIES}_N_/" > $OUTPUT_DIR/${SPECIES}_N.txt

    if [ $SPECIES == "Hcam" ];then
    N="N5"
    elif [ $SPECIES == "Hconf" ];then
    N="N6"        
    elif [ $SPECIES == "Sdef" ];then
    N="N1"
    fi

    col_I=$(head -1 ${hog_dir}/${N}.tsv|tr "\t" "\n"|grep -w -n "${SPECIES}_I"|cut -d ":" -f1)
    col_N=$(head -1 ${hog_dir}/${N}.tsv|tr "\t" "\n"|grep -w -n "${SPECIES}_N"|cut -d ":" -f1)

    echo "SPECIES: $SPECIES $col_I $col_N"
    echo "N: $N"

: <<'END'
>${OUTPUT_DIR}/${SPECIES}_1_1.txt

    while read -r line; do
        gI=$(echo $line|cut -d "|" -f1|sed 's/\s//g'|grep -v '^$'|sed "s/^/${SPECIES}_I_/") 
        gN=$(echo $line|cut -d "|" -f2|sed 's/\s//g'|grep -v '^$'|sed "s/^/${SPECIES}_N_/")

        HOG1=$(grep -w "$gI" ${hog_dir}/${N}.tsv|cut -f1)
        HOG2=$(grep -w "$gN" ${hog_dir}/${N}.tsv|cut -f1)
        # check if HOG1 is empty
        if [ -z "$HOG1" ]; then
            HOG1="NA"
        fi
        # check if HOG2 is empty
        if [ -z "$HOG2" ]; then
            HOG2="NA"
        fi
        echo -e "$gI\t$gN\t$HOG1\t$HOG2" >> $OUTPUT_DIR/${SPECIES}_1_1.txt
    done < $OUTPUT_DIR/${SPECIES}_infected_vs_naive.tsv/gene@mrna@cds_1_1_id_list.txt
END
# D & E: Genes with no genomic overlap between infected and naive AND genes assigned to orthogroup
D=$(grep -o -Fwf $OUTPUT_DIR/${SPECIES}_I.txt ${hog_dir}/${N}.tsv|wc -l)
E=$(grep -o -Fwf $OUTPUT_DIR/${SPECIES}_N.txt ${hog_dir}/${N}.tsv|wc -l)

# A & B: Genes with no genomic overlap between infected and naive AND genes not assigned to any orthogroup (total - D or E)
A=$(( $(wc -l < $OUTPUT_DIR/${SPECIES}_I.txt) - D ))
B=$(( $(wc -l < $OUTPUT_DIR/${SPECIES}_N.txt) - E ))

# C : Genes in both infected and naive AND in same orthogroup ("HOG1" == "HOG2")
C=$(awk -F "\t" '$3==$4' $OUTPUT_DIR/${SPECIES}_1_1.txt|grep -v -w "NA"|wc -l)

# F: Genes in both infected and naive AND not in any orthogroup (both "NA")
F=$(awk -F "\t" '$3=="NA" && $4=="NA"' $OUTPUT_DIR/${SPECIES}_1_1.txt|wc -l)

# Total genes in both infected and naive
T_I=$(( $(wc -l < $OUTPUT_DIR/${SPECIES}_I.txt) + $C + $F ))
T_N=$(( $(wc -l < $OUTPUT_DIR/${SPECIES}_N.txt) + $C + $F ))

echo -e "$SPECIES\t$A\t$B\t$C\t$D\t$E\t$F\t$T_I\t$T_N" >> $OUTPUT_DIR/IvN_gene_summary.txt

# Find out if unique genes are in the immune orthologous gene families



}


export -f get_genes
mkdir -p annot
>annot/IvN_gene_summary.txt
echo -e "SPECIES\tA\tB\tC\tD\tE\tF\tTotal_I\tTotal_N" > annot/IvN_gene_summary.txt
parallel -j 3 --colsep "\t" get_genes {2} :::: species.txt

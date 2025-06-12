#! /bin/bash

while [ "$1" != "" ]; do
    case $1 in
        -DE_genes )            shift
                                DE_genes=$1
                                ;;
        -unique_genes_DEstatus )              shift
                                unique_genes_DEstatus=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done



# This script is used to get the up and down regulated genes from the DESeq2 output file

if [ "$DE_genes" == "true" ]; then
#exit if Orthogroups.tsv and NP_FBid.tsv does not exist
orthofinder_dir="orthofinder/OrthoFinder/Results_Apr05/"
if [ ! -s ${orthofinder_dir}/Orthogroups.tsv ] || [ ! -s ${orthofinder_dir}/NP_FBid.tsv ]; then
    echo "Orthogroups.tsv and NP_FBid.tsv does not exist. Please run orthofinder.sh first"
    exit 1
fi

get_genes () {
    species=$1
    sp=$2
    genome=$3
    orthofinder_dir="orthofinder/OrthoFinder/Results_Apr05/"

    grep -A1 -Fwf <(awk '$3>0' ${sp}_SDE_genes.tsv|cut -f1|sed 's/^/'"${sp}"'_/') <( fasta_formatter -i orthofinder/${sp}.fa )|sed '/^-/d' > DESeq2_genes/${sp}_UP.fasta
    grep -A1 -Fwf <(awk '$3<0' ${sp}_SDE_genes.tsv|cut -f1|sed 's/^/'"${sp}"'_/') <( fasta_formatter -i orthofinder/${sp}.fa )|sed '/^-/d' > DESeq2_genes/${sp}_DOWN.fasta

    # get the Dmel homologs for UP and DOWN genes
    grep ">" DESeq2_genes/${sp}_UP.fasta|cut -d " " -f1|sed 's/>//' > DESeq2_genes/${sp}_DE_genes.txt
    grep ">" DESeq2_genes/${sp}_DOWN.fasta|cut -d " " -f1|sed 's/>//' >> DESeq2_genes/${sp}_DE_genes.txt

    echo -e "Gene\tTranscript\tDmel_Name\tFBgn\tCGid\tHOG\tOG\t$(head -1 ${sp}_SDE_genes.tsv|cut -f2-)" > DESeq2_genes/${sp}_Dmel_genes.tsv

    # get Dmel homologs from the orthofinder.tsv file
    for gene in $(cat DESeq2_genes/${sp}_DE_genes.txt); do
        #remove the sp_ prefix and _t* suffix
        gene_n=$(echo $gene|sed 's/^'"${sp}"'_//;s/\.t[0-9]$//')
        col=$(head -1 ${orthofinder_dir}/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "Dmel" |cut -d ":" -f1)
        sp_col=$(head -1 ${orthofinder_dir}/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${sp}" |cut -d ":" -f1)
        
        HOG=$(cut -f1,$col,$sp_col ${orthofinder_dir}/Orthogroups.tsv|grep -w "$gene"|cut -f1)
        
        if [ -z "$HOG" ]; then
            echo "No HOG found for gene: $gene"
            HOG="NA"
            OG="NA"
            CG="NA"
            FBgn="NA"
            name="NA"
            Dmel_genes=""
        else
            Dmel_genes=$(cut -f1,$col,$sp_col ${orthofinder_dir}/Orthogroups.tsv|grep -w "$HOG"|cut -f2|sed 's/\s//g')
            Dmel_genes=$(echo "$Dmel_genes" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            OG=$(cat ${orthofinder_dir}/Orthogroups.tsv|grep -w "$HOG"|cut -f2)
        fi
        
        if [ -z "$Dmel_genes" ]; then
            echo "No Dmel_genes found for HOG: $gene"
            CG="NA"
            FBgn="NA"
            name="NA"
        else
            CG=$(grep -Fwf <(echo "$Dmel_genes"|tr "," "\n") ${orthofinder_dir}/NP_FBid.tsv|cut -f2|tr "\n" ","|sed 's/,$//')
            FBgn=$(grep -Fwf <(echo "$Dmel_genes"|tr "," "\n")  ${orthofinder_dir}/NP_FBid.tsv|cut -f3|tr "\n" ","|sed 's/,$//')
            name=$(grep -Fwf <(echo "$Dmel_genes"|tr "," "\n")  ${orthofinder_dir}/NP_FBid.tsv|cut -f4|tr "\n" ","|sed 's/,$//')
        fi
        line=$(grep -w "$gene_n" ${sp}_SDE_genes.tsv|cut -f2-)
        echo -e "${gene_n}\t${gene}\t${name}\t${FBgn}\t${CG}\t${HOG}\t${OG}\t${line}" >> DESeq2_genes/${sp}_Dmel_genes.tsv
    done
    
}

export -f get_genes
mkdir -p DESeq2_genes
parallel -j 3 --colsep "\t" get_genes {1} {2} {3} :::: species.txt

fi


if [ "$unique_genes_DEstatus" == "true" ]; then

unique_genes(){
    sp=$1
    orthofinder_dir="orthofinder/OrthoFinder/Results_Apr05/"
    # DE genes that homologs with unique genes
    col_I=$(head -1 ${orthofinder_dir}/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${sp}_I" |cut -d ":" -f1)
    col_N=$(head -1 ${orthofinder_dir}/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${sp}_N" |cut -d ":" -f1)
    grep -Fwf <(grep -Fwf annot/${sp}_I.txt <(grep -Fwf <(cut -f6 DESeq2_genes/${sp}_Dmel_genes.tsv |tail -n +2) ${orthofinder_dir}/Orthogroups.tsv)|awk -F "\t" -v s=$col_N '$s==""'|cut -f1) DESeq2_genes/${sp}_Dmel_genes.tsv |sort -k9,9nr > DESeq2_genes/${sp}_unique_I_genes.txt
    grep -Fwf <(grep -Fwf annot/${sp}_N.txt <(grep -Fwf <(cut -f6 DESeq2_genes/${sp}_Dmel_genes.tsv |tail -n +2) ${orthofinder_dir}/Orthogroups.tsv)|awk -F "\t" -v s=$col_I '$s==""'|cut -f1) DESeq2_genes/${sp}_Dmel_genes.tsv |sort -k9,9nr > DESeq2_genes/${sp}_unique_N_genes.txt

}

export -f unique_genes

parallel -j 3 --colsep "\t" unique_genes {2} :::: species.txt
fi

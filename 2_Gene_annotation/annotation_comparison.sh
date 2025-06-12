#!/bin/bash
# Check if the required arguments are provided
while [ "$1" != "" ]; do
    case $1 in
        -agat )                shift
                                agat=$1
                                ;;
        -table )              shift
                                table=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done
if [ "$agat" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate agat
    echo "Running AGAT for annotation comparison..."

# Run AGAT for each species
run_agat(){
    SPECIES=$1
    GENOME=$2
    OUTPUT_DIR="annot"
    mkdir -p "$OUTPUT_DIR/$SPECIES"

    agat_sp_statistics.pl --gff $OUTPUT_DIR/${SPECIES}/braker.gff3 -o  $OUTPUT_DIR/${SPECIES}/${SPECIES}_stats.tsv -g genomes/"$GENOME" 
    agat_sp_statistics.pl --gff $OUTPUT_DIR/${SPECIES}_infected/braker.gff3 -o $OUTPUT_DIR/${SPECIES}_infected/${SPECIES}_infected_stats.tsv -g genomes/"$GENOME" 
    agat_sp_statistics.pl --gff $OUTPUT_DIR/${SPECIES}_naive/braker.gff3 -o $OUTPUT_DIR/${SPECIES}_naive/${SPECIES}_naive_stats.tsv -g genomes/"$GENOME" 


    # Compare annotations
    agat_sp_compare_two_annotations.pl --gff1 $OUTPUT_DIR/${SPECIES}/braker.gff3 --gff2 $OUTPUT_DIR/${SPECIES}_infected/braker.gff3 --out $OUTPUT_DIR/${SPECIES}_all_vs_infected 
    agat_sp_compare_two_annotations.pl --gff1 $OUTPUT_DIR/${SPECIES}/braker.gff3 --gff2 $OUTPUT_DIR/${SPECIES}_naive/braker.gff3 --out $OUTPUT_DIR/${SPECIES}_all_vs_naive
    agat_sp_compare_two_annotations.pl --gff1 $OUTPUT_DIR/${SPECIES}_infected/braker.gff3 --gff2 $OUTPUT_DIR/${SPECIES}_naive/braker.gff3 --out $OUTPUT_DIR/${SPECIES}_infected_vs_naive

}

export -f run_agat
mkdir -p annot
parallel -j 3 --colsep "\t" run_agat {2} {3} :::: species.txt
fi

# Compare immune genes
if [ "$table" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate agat

    cut -f1,4- orthofinder/OrthoFinder/Results_Apr05/Phylogenetic_Hierarchical_Orthogroups/N0.tsv > orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv
    orthofinder_tsv="orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv"
    #table header
    immuneOG=$(grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) ${orthofinder_tsv}| cut -f1 | sort -u | wc -l)
    #immuneGenes=$(grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) ${orthofinder_tsv}| awk -F'\t' '$2 != ""' |cut -f2 | tr "," "\n" | sort -u | wc -l)
    immuneGenes=$(wc -l immune_rna_Dmel.tsv|cut -d " " -f1)
    # mean length of the immune genes
    grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) ${orthofinder_tsv}| awk -F'\t' '$2 != ""' |cut -f2 | tr "," "\n" | sort -u|sed 's/\s//g' > genes.tmp
    mean_length=$(bioawk -c fastx '{ print $name, length($seq) }' orthofinder/Dmel.fa |grep -Fwf genes.tmp |cut -f2 | ./data_summary.sh|cut -f3)
    
    echo -e "Species\tGene_number\tMean_CDS\tImmune_OGs($immuneOG)\tImmune_genes($immuneGenes)\tmean_length($(echo "$(echo $mean_length) * 0.003" | bc))" > annotation_comparison.tsv

    get_stat(){
        SPECIES=$1
        total_immuneOG=$2
        total_immuneGenes=$3
 
        col=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${SPECIES}"|cut -d ":" -f1)
        # Get gene number and mean CDS length
        stat=$(bioawk -c fastx '{ print $name, length($seq) }' orthofinder/${SPECIES}.fa| cut -f2 | ./data_summary.sh | cut -f2-4)
        GENE_NUMBER=$(echo $stat | cut -d " " -f1)
        MEAN_CDS=$(echo "$(echo $stat | cut -d " " -f2) * 0.003" | bc)

        # Get immune OGs
        IMMUNE_OGS=$(grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != "" {count++} END {print count}')
        grep -Fwf <(grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a == ""'|cut -f2|tr ",\s" "\n"|sed 's/^\s//g') immune_rna_Dmel.tsv > ${SPECIES}_missing_immune_genes.tsv

        grep -Fwf <(grep -Fwvf <(cut -f2 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|tail -n +2|tr "," "\n" |sed 's/\s//g'|sed '/^$/d') <(grep ">" orthofinder/Dmel.fa|sed 's/>//') ) immune_rna_Dmel.tsv >> ${SPECIES}_missing_immune_genes.tsv

        grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a == ""'|cut -f1 > ${SPECIES}_missing_immune_HOGs.tsv
        #n_immune_genes=$(grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != ""' | cut -f${col} | tr "," "\n" | sort -u | wc -l)

        # immune genes= total immune genes - missing immune genes
        n_immune_genes=$(( $total_immuneGenes - $(wc -l ${SPECIES}_missing_immune_genes.tsv| cut -d " " -f1) ))

        # percentage of immune OGs and immune genes
        #PERCENT_IMMUNE_OGS=$(printf "%.2f" "$(echo "scale=4; ($IMMUNE_OGS / $total_immuneOG) * 100" | bc)")
        #PERCENT_IMMUNE_GENES=$(printf "%.2f" "$(echo "scale=4; ($n_immune_genes / $total_immuneGenes) * 100" | bc)")

        # mean length of the immune genes
        grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != ""' | cut -f${col} | tr "," "\n" | sort -u|sed 's/\s//g' > genes.tmp
        mean_length=$(bioawk -c fastx '{ print $name, length($seq) }' orthofinder/${SPECIES}.fa | grep -Fwf genes.tmp | cut -f2 | ./data_summary.sh | cut -f3)
    
        mean_length=$(echo "$(echo $mean_length) * 0.003" | bc)
        
               rm -rf genes.tmp

        echo -e "$SPECIES\t$GENE_NUMBER\t$MEAN_CDS\t$IMMUNE_OGS\t$n_immune_genes\t$mean_length" >> annotation_comparison.tsv

    }

export -f get_stat


for SP in $(cut -f2 species.txt); do
    #set -x

    get_stat ${SP}"_N" $immuneOG $immuneGenes
    get_stat ${SP}"_I" $immuneOG $immuneGenes 
    get_stat ${SP} $immuneOG $immuneGenes


    #set +x
: <<'END'
    # Estimate percentage overlap of genomic coordinates between annotations for immune genes (ALL vs infected vs naive annotations), calculate percentage for each overlap genes (0-100%)
    # Immune gene list (genes.tmp)
    # Get the immune genes from the annotation
    col=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${SP}_I"|cut -d ":" -f1)
    grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != ""' | cut -f${col} | tr "," "\n" | sort -u|sed 's/\s//g' > genes.tmp
    awk '$3=="exon"'  annot/${SP}_infected/braker.gff3|grep -Fwf <(sed "s/${SP}_I_//" genes.tmp) > I_genes.gff

    col=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${SP}_N"|cut -d ":" -f1)
    grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != ""' | cut -f${col} | tr "," "\n" | sort -u|sed 's/\s//g' > genes.tmp
    awk '$3=="exon"'  annot/${SP}_naive/braker.gff3|grep -Fwf <(sed "s/${SP}_N_//" genes.tmp) > N_genes.gff

    col=$(head -1 orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv|tr "\t" "\n"|grep -w -n "${SP}"|cut -d ":" -f1)
    grep -Fwf <(cut -f1 immune_rna_Dmel.tsv) orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv | awk -F'\t' -v a=$col '$a != ""' | cut -f${col} | tr "," "\n" | sort -u|sed 's/\s//g' > genes.tmp
    awk '$3=="exon"'  annot/${SP}/braker.gff3|grep -Fwf <(sed "s/${SP}_//" genes.tmp) > all_genes.gff

    bedtools sort -i I_genes.gff > I_genes.gff.sorted
    bedtools sort -i N_genes.gff > N_genes.gff.sorted
    bedtools sort -i all_genes.gff > all_genes.gff.sorted
    mv I_genes.gff.sorted I_genes.gff
    mv N_genes.gff.sorted N_genes.gff
    mv all_genes.gff.sorted all_genes.gff

    bedtools intersect -a all_genes.gff -b I_genes.gff -names infected -sorted -wao | awk -v comparison="All_vs_infected" 'BEGIN {OFS="\t"} {
        feature_length = $5 - $4 + 1;
        overlap = $NF;
        percent_overlap = (overlap / feature_length) * 100;
        if ($NF == ".") overlap = 0;  # Handle no overlap
        print comparison, $0, percent_overlap
    }' > ${SP}_overlaps.tsv

    bedtools intersect -a all_genes.gff -b N_genes.gff -names naive -sorted -wao | awk -v comparison="All_vs_naive" 'BEGIN {OFS="\t"} {
        feature_length = $5 - $4 + 1;
        overlap = $NF;
        percent_overlap = (overlap / feature_length) * 100;
        if ($NF == ".") overlap = 0;  # Handle no overlap
        print comparison, $0, percent_overlap
    }' >> ${SP}_overlaps.tsv

    bedtools intersect -a I_genes.gff -b N_genes.gff -names naive -sorted -wao | awk -v comparison="infected_vs_naive" 'BEGIN {OFS="\t"} {
        feature_length = $5 - $4 + 1;
        overlap = $NF;
        percent_overlap = (overlap / feature_length) * 100;
        if ($NF == ".") overlap = 0;  # Handle no overlap
        print comparison, $0, percent_overlap
    }' >> ${SP}_overlaps.tsv

    rm -rf genes.tmp I_genes.gff N_genes.gff all_genes.gff
    echo "$SP  done"

END

done


# get the table of immune genes present/absent in all species
grep -a -Fwf <(cut -f3 immune_rna_Dmel.tsv) Immune_genes_categories_MarkHanson.txt > immune.tmp

echo -e "Hcam\tHconf\tSdef\tHOG\t$(head -1 Immune_genes_categories_MarkHanson.txt)" > immune_genes_summary.tsv
echo "start"

while IFS=$'\t' read -r line; do
    # Find out if the gene is missing in the species
    gene=$(echo "$line" | awk -F'\t' '{print $3}')
    
    Hcam_status="Yes"
    Hconf_status="Yes"
    Sdef_status="Yes"

    rna=$(grep -w "$gene" immune_rna_Dmel.tsv |awk -F'\t' '{print $1}')
    if [ -z "$rna" ]; then
        HOG="NA"
    else
        HOG=$(grep -w "$rna" orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv |awk -F'\t' '{print $1}')
        # 
        if [ -z "$HOG" ]; then
            HOG="NA"
        fi
    fi

    Hcam=$(grep -w "$gene" <(grep -Fwvf <(cut -f9 missing_genes/Hcam_30_tophits.gff |cut -d ";" -f2|sed -e 's/sequence//' -e 's/\s//g') Hcam_missing_immune_genes.tsv))
    Hconf=$(grep -w "$gene" <(grep -Fwvf <(cut -f9 missing_genes/Hconf_30_tophits.gff |cut -d ";" -f2|sed -e 's/sequence//' -e 's/\s//g') Hconf_missing_immune_genes.tsv))
    Sdef=$(grep -w "$gene" <(grep -Fwvf <(cut -f9 missing_genes/Sdef_30_tophits.gff |cut -d ";" -f2|sed -e 's/sequence//' -e 's/\s//g') Sdef_missing_immune_genes.tsv))

    if [ -n "$Hcam" ]; then
        Hcam_status="No"
    fi
    if [ -n "$Hconf" ]; then
        Hconf_status="No"
    fi
    if [ -n "$Sdef" ]; then
        Sdef_status="No"
    fi

    echo -e "$Hcam_status\t$Hconf_status\t$Sdef_status\t$HOG\t$line" >> immune_genes_summary.tsv
done < immune.tmp

rm -rf immune.tmp


fi
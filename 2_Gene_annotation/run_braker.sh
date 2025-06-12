#!/bin/bash

# Check if the required arguments are provided
while [ "$1" != "" ]; do
    case $1 in
        -exonerate )            shift
                                exonerate=$1
                                ;;
        -braker )              shift
                                braker=$1
                                ;;
        -braker2 )             shift
                                braker2=$1
                                ;;
        -addUTR )              shift
                                addUTR=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done


if [ "$exonerate" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate annotation

    echo "Running Exonerate for hint generation..."
#run exonerate to generate hints for BRAKER3
# Split Dmel.fa into ~100 smaller files evenly
split_fasta_evenly() {
    INPUT_FASTA=$1
    OUTPUT_DIR=$2
    CHUNKS=$3

    mkdir -p "$OUTPUT_DIR"
    seqkit split2 -p "$CHUNKS" -O "$OUTPUT_DIR" "$INPUT_FASTA"
}
export -f split_fasta_evenly

# Run Exonerate on each split file and genome
run_exonerate_parallel() {
    SPLIT_DIR=$1
    GENOME=$2
    OUTPUT_DIR=$3
    SPECIES=$4

    # Run Exonerate in parallel for each split file
    find "$SPLIT_DIR" -name "*.fa" | parallel -j 40 exonerate --model protein2genome --percent 30 {} genomes/${GENOME} \
        --showalignment FALSE --showvulgar FALSE --showtargetgff yes '>' ${OUTPUT_DIR}/${SPECIES}/{/.}_${SPECIES}_exonerate.gff
}


    # Main script
PROTEINS="Dmel.fa"
GENOME_DIR="genomes"
OUTPUT_DIR="exonerate_hints"
SPLIT_DIR="Dmel_proteins"
CHUNKS=100  # Number of files to split into

# Step 1: Split the Dmel.fa file evenly
# export -f split_fasta_evenly
#split_fasta_evenly "$PROTEINS" "$SPLIT_DIR" "$CHUNKS"

# Step 2: Run Exonerate for each species genome
while read -r FullName SPECIES GENOME; do
    export -f run_exonerate_parallel
    mkdir -p "$OUTPUT_DIR/$SPECIES"
    run_exonerate_parallel "$SPLIT_DIR" "$GENOME" "$OUTPUT_DIR" "$SPECIES"
    # merge all the gff files into one
    cat ${OUTPUT_DIR}/${SPECIES}/*_${SPECIES}_exonerate.gff > ${OUTPUT_DIR}/${SPECIES}_exonerate.gff
    
done < species.txt

    # Convert Exonerate output to hints
    #exonerate2hints.pl --in ${OUTPUT_DIR}/${sp}_exonerate.gff --minintronlen=10 --CDSpart_cutoff=5 --out=${OUTPUT_DIR}/${sp}_hints.gff

fi


if [ "$braker" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate braker3

    echo "Running BRAKER3 for gene prediction..."

    # Run BRAKER3 for each species (infected vs naive RNAseq)
    run_braker(){
        SPECIES=$1
        GENOME=$2
        OUTPUT_DIR="annot"

        mkdir -p "$OUTPUT_DIR/${SPECIES}_infected"
	rm -rf ~/fly_annotation/tools_braker3/Augustus/bin/../config/species/"$SPECIES".infected
        #braker.pl --species="$SPECIES".infected --genome=genomes/"$GENOME" --prot_seq="Dmel.fa"  --bam=mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam --threads=20 --workingdir="${OUTPUT_DIR}/${SPECIES}_infected" --gff3 --softmasking --AUGUSTUS_ab_initio --augustus_args="--species=fly"
        gushr.py -t ${OUTPUT_DIR}/${SPECIES}_infected/braker.gtf -b mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam -g genomes/"$GENOME" -o ${SPECIES}_infected_utrs -c 20
        mv ${SPECIES}_infected_utrs* ${OUTPUT_DIR}/${SPECIES}_infected/
	echo "$SPECIES infected done"

	mkdir -p "$OUTPUT_DIR/${SPECIES}_naive"
	rm -rf ~/fly_annotation/tools_braker3/Augustus/bin/../config/species/"$SPECIES".naive
	#braker.pl --species="$SPECIES".naive --genome=genomes/"$GENOME" --prot_seq="Dmel.fa"  --bam=mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam --threads=20 --workingdir="${OUTPUT_DIR}/${SPECIES}_naive" --gff3 --softmasking --AUGUSTUS_ab_initio --augustus_args="--species=fly"
    gushr.py -t ${OUTPUT_DIR}/${SPECIES}_naive/braker.gtf -b mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam -g genomes/"$GENOME" -o ${SPECIES}_naive_utrs -c 20
	mv ${SPECIES}_naive_utrs* ${OUTPUT_DIR}/${SPECIES}_naive/
    echo "$SPECIES naive done"
    }

    export -f run_braker
    mkdir -p annot
    parallel -j 3 --colsep "\t" run_braker {2} {3} :::: species.txt

fi


if [ "$braker2" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate braker3

    echo "Running BRAKER3 for gene prediction..."
    # Run BRAKER3 for each species
    run_braker(){
        SPECIES=$1
        GENOME=$2
        OUTPUT_DIR="annot"
        mkdir -p "$OUTPUT_DIR/$SPECIES"

        #braker.pl --species="$SPECIES".1 --genome=genomes/"$GENOME" --prot_seq="Dmel.fa"  --bam=mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam,mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam --threads=24 --workingdir="${OUTPUT_DIR}/${SPECIES}" --gff3 --softmasking --AUGUSTUS_ab_initio --augustus_args="--species=fly"
        # add UTRs
        samtools merge mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam -o mapped_reads/${SPECIES}_merged.bam --threads 24
        gushr.py -t ${OUTPUT_DIR}/${SPECIES}/braker.gtf -b mapped_reads/${SPECIES}_merged.bam -g genomes/"$GENOME" -o ${SPECIES}_utrs -c 24
        mv ${SPECIES}_utrs* ${OUTPUT_DIR}/${SPECIES}/
        rm -rf mapped_reads/${SPECIES}_merged.bam
    }
    export -f run_braker
    mkdir -p annot
    parallel -j 3 --colsep "\t" run_braker {2} {3} :::: species.txt

fi


if [ "$addUTR" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate braker3

    add_UTRs(){
        SPECIES=$1
        GENOME=$2
        OUTPUT_DIR="annot"
        mkdir -p "$OUTPUT_DIR/$SPECIES"

        gushr.py -t ${OUTPUT_DIR}/${SPECIES}_infected/braker.gtf -b mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam -g genomes/"$GENOME" -o ${SPECIES}_infected_utrs -c 40
        mv ${SPECIES}_infected_utrs* ${OUTPUT_DIR}/${SPECIES}_infected/
	    echo "$SPECIES infected done"
	
        gushr.py -t ${OUTPUT_DIR}/${SPECIES}_naive/braker.gtf -b mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam -g genomes/"$GENOME" -o ${SPECIES}_naive_utrs -c 40
	    mv ${SPECIES}_naive_utrs* ${OUTPUT_DIR}/${SPECIES}_naive/
        echo "$SPECIES naive done"
        # add UTRs
        samtools merge mapped_reads/${SPECIES}_infected_Aligned.sortedByCoord.out.bam mapped_reads/${SPECIES}_naive_Aligned.sortedByCoord.out.bam -o mapped_reads/${SPECIES}_merged.bam --threads 40
        gushr.py -t ${OUTPUT_DIR}/${SPECIES}/braker.gtf -b mapped_reads/${SPECIES}_merged.bam -g genomes/"$GENOME" -o ${SPECIES}_utrs -c 40
        mv ${SPECIES}_utrs* ${OUTPUT_DIR}/${SPECIES}/
        rm -rf mapped_reads/${SPECIES}_merged.bam
    }

    export -f add_UTRs
    mkdir -p annot
    parallel -j 1 --colsep "\t" add_UTRs {2} {3} :::: species.txt

fi


#nohup orthofinder -f orthofinder/ -t 60 -s annot/SpeciesTree_rooted.txt -a 12 -M msa -S blast -T iqtree -A mafft -X  > orthofinder.nohup.log 2>&1 &

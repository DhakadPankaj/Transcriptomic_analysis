#! /bin/bash

# This script is used to run the de novo assembly of unmapped reads
# using rnaSPAdes.

# Extract unmapped reads from the bam files

while [ "$1" != "" ]; do
    case $1 in
        -UnmappedReads )            shift
                                UnmappedReads=$1
                                ;;
        -deNovoAssembly )              shift
                                deNovoAssembly=$1
                                ;;
        -TaxonomicClassification )        shift
                                TaxonomicClassification=$1
                                ;;
        -Remapping )              shift
                                Remapping=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done

if [ "$UnmappedReads" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate annotation
    mkdir -p unmapped_reads

get_unmapped () {
    sp=$1
    treat=$2
    lable=$3
    threads=10
    mkdir -p trimmed_reads/${treat}/${sp}
    #zcat trimmed_reads/${treat}/${sp}/${lable}_*_1.fq.gz > trimmed_reads/${treat}/${sp}/${lable}_1.fastq
    #zcat trimmed_reads/${treat}/${sp}/${lable}_*_2.fq.gz > trimmed_reads/${treat}/${sp}/${lable}_2.fastq

    #STAR --runThreadN $threads --genomeDir genomes/${sp}_STAR_gff --readFilesIn trimmed_reads/${treat}/${sp}/${lable}_1.fastq trimmed_reads/${treat}/${sp}/${lable}_2.fastq --outFileNamePrefix unmapped_reads/bam/${sp}_${lable}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outSAMunmapped Within

    #samtools index unmapped_reads/bam/${sp}_${lable}_Aligned.sortedByCoord.out.bam

    #gunzip unmapped_reads/bam/${sp}_${lable}_Aligned.sortedByCoord.out.bam.gz
    samtools view -b -f 12 -F 256 unmapped_reads/bam/${sp}_${lable}_Aligned.sortedByCoord.out.bam > unmapped_reads/${sp}_${lable}_unmapped.bam
    
    samtools sort -n -@ $threads unmapped_reads/${sp}_${lable}_unmapped.bam -o unmapped_reads/${sp}_${lable}_unmapped_sorted.bam

    samtools fastq -1 unmapped_reads/${sp}_${lable}_unmapped_R1.fastq.gz -2 unmapped_reads/${sp}_${lable}_unmapped_R2.fastq.gz -0 /dev/null -s /dev/null -n unmapped_reads/${sp}_${lable}_unmapped_sorted.bam --threads $threads
    
    rm -rf unmapped_reads/${sp}_${lable}_unmapped.bam
    #gzip unmapped_reads/bam/${sp}_${lable}_Aligned.sortedByCoord.out.bam
    #rm -rf trimmed_reads/${treat}/${sp}/${lable}_1.fastq trimmed_reads/${treat}/${sp}/${lable}_2.fastq
}

export -f get_unmapped
mkdir -p unmapped_reads/bam
# genome index for STAR
#STAR --runThreadN 20 --genomeDir "genomes/${sp}_STAR_gff/" --runMode genomeGenerate --sjdbGTFfile annot/${sp}/braker.gtf --limitGenomeGenerateRAM 150000000000 --genomeFastaFiles genomes/${genome} --genomeSAindexNbases 12 --sjdbOverhang 147
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt |cut -f2,5,7 |parallel -j 5 --colsep "\t" get_unmapped {1} {2} {3}

fi


if [ "$deNovoAssembly" == "true" ]; then

source ~/miniconda3/etc/profile.d/conda.sh
conda activate spades

run_spade(){
    sp=$1
    treat=$2
    lable=$3
    threads=12
    
    mkdir -p unmapped_reads/deNovo_assembly/${sp}
    cd unmapped_reads/deNovo_assembly/${sp}
    rnaspades.py -1 ../../${sp}_${lable}_unmapped_R1.fastq.gz -2 ../../${sp}_${lable}_unmapped_R2.fastq.gz -o ${lable} --threads $threads
    cd -
}
export -f run_spade
mkdir -p unmapped_reads/deNovo_assembly/
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | cut -f2,5,7 |parallel -j 5 --colsep "\t" run_spade {1} {2} {3}
fi

# diamond blastp to get taxonomic classification
if [ "$TaxonomicClassification" == "true" ]; then

run_diamond(){
threads=15
sp=$1
treat=$2
lable=$3

: <<'END'
# getORF
cd unmapped_reads/deNovo_assembly/${sp}
getorf -sequence ${lable}/transcripts.fasta -minsize 600  -outseq ${lable}/orfs.fasta

# diamond blastp
diamond blastp -d /data/BLAST_databases/nr_diamond -q ${lable}/orfs.fasta -o ${lable}/diamond_results.tsv --max-target-seqs 5 --evalue 1e-20 --outfmt 6 qseqid sseqid pident length evalue bitscore stitle sskingdoms skingdoms --threads $threads --include-lineage

cd -

END

#prepend lable to the output file
#awk -FS"\t"  -v sample="${lable}" -v sp="${sp}_${treat}" '{OFS="\t"; $1=sample"_"$1; print sp, $0}' unmapped_reads/deNovo_assembly/${sp}/${lable}/diamond_results.tsv | awk -F "\t" -v OFS="\t" '{ if ($8 ~ /\[.*\]/) { match($8, /\[([^\]]+)\]/, arr); print $0, arr[1] } }'|awk -F "\t" -v OFS="\t" '!seen[$2]++' > unmapped_reads/diamond_results/${sp}_${lable}_hits.tsv
awk -F "\t" -v sample="${lable}" -v sp="${sp}_${treat}" '{
    OFS="\t"; 
    $1 = sample"_"$1; 
    print sp, $0
}' unmapped_reads/deNovo_assembly/${sp}/${lable}/diamond_results.tsv | \
awk -F "\t" -v OFS="\t" '{
    if ($8 ~ /\[.*\]/) { 
        match($8, /\[([^\]]+)\]$/, arr);  # Match the last set of square brackets
        print $0, arr[1];  # Append the extracted species name
    } 
}' | \
awk -F "\t" -v OFS="\t" '!seen[$2]++' > unmapped_reads/diamond_results/${sp}_${lable}_hits.tsv


}

export -f run_diamond
mkdir -p unmapped_reads/diamond_results
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | cut -f2,5,7 |parallel -j 5 --colsep "\t" run_diamond {1} {2} {3}

cat unmapped_reads/diamond_results/*_hits.tsv > unmapped_reads/All_diamond_results.tsv

awk -F "\t" '{print $11}' unmapped_reads/All_diamond_results.tsv |sort -u|grep -w -v -e "Bacteria$" -e "Serratia$" -e "Yersinia$" >taxa_list.txt

taxonkit name2taxid taxa_list.txt > name2taxid.txt

cut -f2 name2taxid.txt |taxonkit lineage |taxonkit reformat -r Unassigned -F --format "{p};{c};{o};{f};{g};{s}"|awk -F "\t" -v OFS="\t" '
{
    if ($2 ~ /^cellular organisms;/) {
        split($2, arr, ";");  # Split column 2 by semicolon
        $3 = arr[2] ";" $3;  # Add the second value (e.g., Eukaryota) to column 3
    } else if ($2 ~ /^Viruses;/) {
        $3 = "Viruses;" $3;  # Add "Viruses;" to column 3
    }  else {$3= "unclassified;" $3;}
    print $0;  # Print the modified line
}'|cut -d $'\t' -f1,3 > taxa_classification.txt


fi


# Remapping the reads to the rna assembly
# This is done to get the relative abundance of genus (metagenomics) in the unmapped reads
if [ "$Remapping" == "true" ]; then
source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

run_remap(){
    sp=$1
    treat=$2
    lable=$3
    threads=15
    mkdir -p unmapped_reads/remapping/${sp}
: <<'END' 
    # genome index for STAR
    genome="unmapped_reads/deNovo_assembly/${sp}/${lable}"

    rm -rf ${genome}/tmp/
    STAR --runThreadN $threads --genomeDir "${genome}/${label}_STAR/" --outTmpDir "${genome}/tmp/"  --runMode genomeGenerate --limitGenomeGenerateRAM 150000000000 --genomeFastaFiles ${genome}/transcripts.fasta --genomeSAindexNbases 10

    zcat unmapped_reads/${sp}_${lable}_unmapped_R1.fastq.gz > unmapped_reads/${sp}_${lable}_unmapped_R1.fastq
    zcat unmapped_reads/${sp}_${lable}_unmapped_R2.fastq.gz > unmapped_reads/${sp}_${lable}_unmapped_R2.fastq
    
    STAR --runThreadN $threads --genomeDir "${genome}/${label}_STAR/" --readFilesIn unmapped_reads/${sp}_${lable}_unmapped_R1.fastq unmapped_reads/${sp}_${lable}_unmapped_R2.fastq --outFileNamePrefix unmapped_reads/remapping/${sp}/${lable}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outSAMunmapped Within

    samtools index unmapped_reads/remapping/${sp}/${lable}_Aligned.sortedByCoord.out.bam
    samtools idxstats unmapped_reads/remapping/${sp}/${lable}_Aligned.sortedByCoord.out.bam > unmapped_reads/remapping/${sp}/${lable}_contigs.idxstats

    # $label at the start of contig name
    awk -F "\t" -v sample="${lable}" -v sp="${sp}_${treat}" '{
    OFS="\t";
    $1 = sample"_"$1;
    print sp, $0
    }' unmapped_reads/remapping/${sp}/${lable}_contigs.idxstats > unmapped_reads/remapping/${sp}_${lable}_contigs.idxstats.tmp


    total_reads=$(($(cat unmapped_reads/${sp}_${lable}_unmapped_R1.fastq |wc -l)/4))

    echo -e "${sp}_${treat}\t${lable}\t${total_reads}" > unmapped_reads/remapping/${sp}_${lable}_total_reads.tmp

    rm -rf unmapped_reads/${sp}_${lable}_unmapped_R1.fastq unmapped_reads/${sp}_${lable}_unmapped_R2.fastq
END

# totoal reads mapped to fly genome
reads=$(samtools idxstats unmapped_reads/bam/${sp}_${lable}_Aligned.sortedByCoord.out.bam |cut -f3|./data_summary.sh |cut -f1)

echo -e "${lable}\t${reads}" > unmapped_reads/bam/${sp}_${lable}_total_reads.tmp

}


export -f run_remap
mkdir -p unmapped_reads/remapping
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | cut -f2,5,7 |parallel -j 5 --colsep "\t" run_remap {1} {2} {3}

cat unmapped_reads/remapping/*_contigs.idxstats.tmp > unmapped_reads/All_contigs.idxstats.tsv
cat unmapped_reads/remapping/*_total_reads.tmp > unmapped_reads/All_total_reads.tsv
cat unmapped_reads/bam/*_total_reads.tmp > unmapped_reads/fly_mapped_read_counts.tsv
fi

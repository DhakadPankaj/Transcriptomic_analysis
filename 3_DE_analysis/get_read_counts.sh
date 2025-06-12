while [ "$1" != "" ]; do
    case $1 in
        -annotation )            shift
                                annotation=$1
                                ;;
        -mapping )              shift
                                mapping=$1
                                ;;
        -featureCounts )        shift
                                featureCounts=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done


if [ "$annotation" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate braker3
    run_annot () {
        species=$1
        sp=$2
        genome=$3
        mkdir -p annot/${sp}
        braker.pl --species=${sp} --genome=genomes/${genome} --AUGUSTUS_ab_initio --threads=20 --gff3 --workingdir=annot/${sp} --augustus_args="--species=fly" --bam=mapped_reads/${sp}_infected_Aligned.sortedByCoord.out.bam,mapped_reads/${sp}_naive_Aligned.sortedByCoord.out.bam
    }
    
    export -f run_annot
    mkdir -p annot
    parallel -j 3 --colsep "\t" -d "\r\n"  run_annot {1} {2} {3} :::: species.txt
fi


run_mapping () {
    sp=$1
    treat=$2
    lable=$3
    threads=8
    
    mkdir -p trimmed_reads/${treat}/${sp}
    zcat trimmed_reads/${treat}/${sp}/${lable}_*_1.fq.gz > trimmed_reads/${treat}/${sp}/${lable}_1.fastq
    zcat trimmed_reads/${treat}/${sp}/${lable}_*_2.fq.gz > trimmed_reads/${treat}/${sp}/${lable}_2.fastq

    STAR --runThreadN $threads --genomeDir genomes/tmp/${sp}_STAR_gff --readFilesIn trimmed_reads/${treat}/${sp}/${lable}_1.fastq trimmed_reads/${treat}/${sp}/${lable}_2.fastq --outFileNamePrefix mapped_reads/DFE2/${sp}_${lable}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --quantMode TranscriptomeSAM GeneCounts

    samtools index mapped_reads/DFE2/${sp}_${lable}_Aligned.sortedByCoord.out.bam
    
    rm -rf trimmed_reads/${treat}/${sp}/${lable}_1.fastq trimmed_reads/${treat}/${sp}/${lable}_2.fastq
}


# only run this if -mapping is specified
if [ "$mapping" == "true" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate annotation
    export -f run_mapping
    mkdir -p mapped_reads/DFE2
 mkdir -p genomes/tmp/
  while IFS=$'\t' read species sp genome; do

  # genome index for STAR
  STAR --runThreadN 20 --genomeDir "genomes/tmp/${sp}_STAR_gff/" --runMode genomeGenerate --sjdbGTFfile annot/${sp}/braker.gtf --limitGenomeGenerateRAM 150000000000 --genomeFastaFiles genomes/${genome} --genomeSAindexNbases 12 --sjdbOverhang 147

  tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="infected"'|cut -f2,5,7 |parallel -j 5 --colsep "\t" run_mapping {1} {2} {3}
  tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="naive"'|cut -f2,5,7 |parallel -j 5  --colsep "\t" run_mapping {1} {2} {3}
  done < tmp.txt
  
fi



#  if featureCounts is true

if [ "$featureCounts" == "true" ]; then
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
    run_featureCounts () {
        sp=$1
        treat=$2
        lable=$3
        featureCounts -T 16 -p -t gene -g ID -a annot/${sp}/braker.gff3 -o counts/${lable}.gene.txt mapped_reads/DFE/${sp}_${lable}_Aligned.sortedByCoord.out.bam
    }
    
    export -f run_featureCounts
    mkdir -p counts
    while IFS=$'\t' read species sp genome; do
        tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="infected"'|cut -f2,5,7 |parallel -j 5 --colsep "\t" run_featureCounts {1} {2} {3}
        tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="naive"'|cut -f2,5,7 |parallel -j 5  --colsep "\t" run_featureCounts {1} {2} {3}

        labels=$(tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a'|cut -f7|sort|uniq)

        paste_cmd="paste"
        for label in $labels; do
            # for the first label keep Geneid
            if [ "$paste_cmd" == "paste" ]; then
                paste_cmd+=" <(awk -v lable=${label} 'BEGIN {OFS=\"\t\"} NR==1 {print \"Geneid\", lable} NR>2 {print \$1, \$7}' counts/${label}.gene.txt)"
                else
                paste_cmd+=" <(awk -v lable=${label} 'BEGIN {OFS=\"\t\"} NR==1 {print lable} NR>2 {print \$7}' counts/${label}.gene.txt)"
            fi
        done
        eval $paste_cmd > counts/${sp}_combined_counts.txt

    done < species.txt



fi

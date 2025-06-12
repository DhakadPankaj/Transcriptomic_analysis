#! /bin/bash

# get arguments
# Usage: bash run_fastp.sh -trimming -mapping 

source ~/miniconda3/etc/profile.d/conda.sh
conda activate annotation

while [ "$1" != "" ]; do
    case $1 in
        -trimming )            shift
                                trimming=$1
                                ;;
        -mapping )              shift
                                mapping=$1
                                ;;
        -annotation )           shift
                                annotation=$1
                                ;;
        * )                     echo "Invalid argument"
                                exit 1
    esac
    shift
done

# Run fastp to trim reads
run_trim () {
    lable=$1
    lane=$2

    #fastqc -o fastqc_output RawData/${lable}/${lable}_${lane}_1.fq.gz
    #fastqc -o fastqc_output RawData/${lable}/${lable}_${lane}_2.fq.gz
    accession=${lable}_${lane}

    if [ ! -s trimmed_reads/${accession}_1.fq.gz ]; then
      echo "Trimming ${accession}"
      fastp -i RawData/${lable}/${lable}_${lane}_1.fq.gz \
      -I RawData/${lable}/${lable}_${lane}_2.fq.gz \
      -o trimmed_reads/${accession}_1.fq.gz \
      -O trimmed_reads/${accession}_2.fq.gz \
      -h qc/${accession}.html \
      -j qc/${accession}.json \
      -c -y -Y 20 \
      --thread 8
    fi
    
}


# only run this if -trimming is specified
if [ "$trimming" == "true" ]; then
  mkdir -p qc
  mkdir -p trimmed_reads
  export -f run_trim
  tail -n +2 DataQualitySummary.txt | parallel -j 8 --colsep "\t" run_trim {1} {2}
fi


# Run STAR to map reads

run_mv () {
    sp=$1
    treat=$2
    lable=$3
    
    mkdir -p trimmed_reads/${treat}/${sp}
    mv trimmed_reads/${lable}_*.fq.gz trimmed_reads/${treat}/${sp}/

}

# only run this if -mapping is specified
if [ "$mapping" == "true" ]; then
  mkdir -p mapped_reads
  ## comment
  echo "I'm here"

  threads=40
  mem=300G

  export -f run_mv
  while IFS=$'\t' read species sp genome; do
    echo -e "Mapping ${sp} ${genome}"
    mkdir -p trimmed_reads/infected
    mkdir -p trimmed_reads/naive
    mkdir -p trimmed_reads/norm
: <<'END'
    
    tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="infected"'|cut -f2,5,7 |parallel -j 4 --colsep "\t" run_mv {1} {2} {3}
    tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="naive"'|cut -f2,5,7 |parallel -j 4  --colsep "\t" run_mv {1} {2} {3}

    
    # Generate genome index for STAR
    #STAR --runThreadN $threads --genomeDir "genomes/${sp}_STAR/" --runMode genomeGenerate --limitGenomeGenerateRAM 150000000000 --genomeFastaFiles genomes/${genome} --genomeSAindexNbases 12
END

    for treat in "infected" "naive" ; do

#: <<'END2'

#    zcat trimmed_reads/${treat}/${sp}/*_1.fq.gz > trimmed_reads/${treat}/${sp}/${sp}_${treat}_R1.fastq
#    zcat trimmed_reads/${treat}/${sp}/*_2.fq.gz > trimmed_reads/${treat}/${sp}/${sp}_${treat}_R2.fastq
    
   bbnorm.sh in="trimmed_reads/${treat}/${sp}/${sp}_${treat}_R1.fastq" in2="trimmed_reads/${treat}/${sp}/${sp}_${treat}_R2.fastq" out="trimmed_reads/norm/${sp}_${treat}_norm_R1.fastq" out2="trimmed_reads/norm/${sp}_${treat}_norm_R2.fastq" ignorebadquality qin=auto qout=auto min=2 target=300 prefilter threads="$threads" -Xmx"$mem"

    #accession=${lable}_${lane}
    STAR --runThreadN $threads --genomeDir genomes/${sp}_STAR --readFilesIn trimmed_reads/norm/${sp}_${treat}_norm_R1.fastq trimmed_reads/norm/${sp}_${treat}_norm_R2.fastq --outSAMstrandField intronMotif --outFileNamePrefix mapped_reads/${sp}_${treat}_ --outSAMtype BAM SortedByCoordinate 
    #samtools index 

    done
#END
  done < species.txt

fi

if [ "$annotation" == "true" ]; then
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate braker3

  run_braker () {
    species=$1
    sp=$2
    genome=$3

    echo -e "Annotating ${sp} ${genome}"

    for treat in "infected" "naive" ; do
      # Run Braker
    samtools index mapped_reads/${sp}_${treat}_Aligned.sortedByCoord.out.bam
    braker.pl --species=fly --genome=genomes/${genome} --AUGUSTUS_ab_initio --threads=16 --useexisting --gff3 --workingdir=annot/${sp}_${treat} --augustus_args="--species=fly" --bam=mapped_reads/${sp}_${treat}_Aligned.sortedByCoord.out.bam
    done

    # both infected and naive
    #braker.pl --species=fly --genome=genomes/${genome} --AUGUSTUS_ab_initio --threads=20 --useexisting --gff3 --workingdir=annot/${sp} --augustus_args="--species=fly" --bam=mapped_reads/${sp}_Aligned.sortedByCoord.out.bam

  }

  export -f run_braker
  mkdir -p annot
  parallel -j 3 --colsep "\t" -d "\r\n"  run_braker {1} {2} {3} :::: species.txt

fi

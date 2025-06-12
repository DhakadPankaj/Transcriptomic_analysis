#! /bin/bash

# This script will generate a summary of the alignment statistics for the fly infection data



# Function to get the alignment summary
get_summary(){
    # Directory for the fly infection mapped data
    data_dir="/data/home/s2215768/fly_infection_data/mapped_reads/DFE"
    # Output file for the alignment summary
    output_file="/data/home/s2215768/fly_infection_data/alignment_summary.tsv"
    sp=$1
    treat=$2
    lable=$3
    #get stats from Log.final.out
    total_reads=$(grep "Number of input reads" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    mapped_reads=$(grep "Uniquely mapped reads number" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    percent_mapped=$(grep "Uniquely mapped reads %" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    p_multiple_mapped=$(grep "% of reads mapped to multiple loci" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    p_toomany_loci=$(grep "% of reads mapped to too many loci" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    p_toomany_mismatch=$(grep "% of reads unmapped: too many mismatches" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    p_tooshort=$(grep "% of reads unmapped: too short" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)
    p_other=$(grep "% of reads unmapped: other" ${data_dir}/${sp}_${lable}_Log.final.out | cut -f2)

    # Print the alignment summary
    echo -e "$sp\t$treat\t$lable\t$total_reads\t$mapped_reads\t$percent_mapped\t$p_multiple_mapped\t$p_toomany_loci\t$p_toomany_mismatch\t$p_tooshort\t$p_other" >> $output_file
}

# Create a header for the output file
echo -e "species\ttreatment\tlabel\ttotal_reads\tmapped_reads\tpercent_mapped\tp_multiple_mapped\tp_toomany_loci\tp_toomany_mismatch\tp_tooshort\tp_other" > /data/home/s2215768/fly_infection_data/alignment_summary.tsv
export -f get_summary
# Loop over the species
for sp in "Hconf" "Hcam" "Sdef"
do 
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="infected"'|cut -f2,5,7 |parallel -j 1 --colsep "\t" get_summary {1} {2} {3}
tail -n +2 Fly_infection_Prettgeri_WithSubmissionNames.txt | awk -F "\t" -v a=${sp} '$2==a && $5=="naive"'|cut -f2,5,7 |parallel -j 1  --colsep "\t" get_summary {1} {2} {3}
done
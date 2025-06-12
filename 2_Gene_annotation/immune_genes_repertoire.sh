

blastp -query All_Dmel_immune_genes.fa -db orthofinder/Sdef.fa -out Sdef.blastp -evalue 1e-20 -num_threads 16 -outfmt 6 -max_target_seqs 1
blastp -query All_Dmel_immune_genes.fa -db orthofinder/Hconf.fa -out Hconf.blastp -evalue 1e-20 -num_threads 16 -outfmt 6 -max_target_seqs 1
blastp -query All_Dmel_immune_genes.fa -db orthofinder/Hcam.fa -out Hconf.blastp -evalue 1e-20 -num_threads 16 -outfmt 6 -max_target_seqs 1

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
        
        if [ -z "$HOG" ]; then
            HOG="NA"
        else
            sp_genes=$(grep -w "$HOG" orthofinder/OrthoFinder/Results_Apr05/Orthogroups.tsv |awk -F'\t' '{print $2}')
            
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
#! /bin/bash

# Run eggNOG-mapper for each species
run_emapper(){
    SPECIES=$1
    OUTPUT_DIR="eggNOG"
    mkdir -p "$OUTPUT_DIR/$SPECIES"
    emapper.py -o ${SPECIES} -m diamond --cpu 30 --output_dir $OUTPUT_DIR/$SPECIES --override -m hmmer -d Diptera -i orthofinder/${SPECIES}.fa --evalue 0.05 --score 60 --pident 40 --query_cover 20 --subject_cover 30 --itype proteins --tax_scope 7147 --target_orthologs all --go_evidence non-electronic --pfam_realign realign --report_orthologs --decorate_gff yes --excel
}
    
export -f run_emapper
mkdir -p eggNOG
parallel -j 3 --colsep "\t" run_emapper {2} :::: species.txt
parallel -j 3 --colsep "\t" run_emapper {2}"_I" :::: species.txt
parallel -j 3 --colsep "\t" run_emapper {2}"_N" :::: species.txt
}

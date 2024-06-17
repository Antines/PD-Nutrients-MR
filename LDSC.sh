## Linkage disequilibrium score regression (LDSR) between nutrients and PD

#All GWAS full summary statistics files universally processed
    awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print "SNP", "A1", "A2", "Zscore", "Pvalue", "N"; next} {print $<SNPcol>, $<Effect Allele>, $<Other Allele>, $<Beta>/$<SE>, $<Pvalue, $<Sample size>}' summary_stats_file.tsv > summary_stats_ldsc.tsv

#!/bin/bash
mkdir -p Munged
mkdir -p PD_Sup_results

for FILE in ./summary_stats/*ldsc.tsv; do
    BASENAME=$(basename $FILE)
    # Run munge sumstats
    echo "Munging summary statistics for $BASENAME"
    ../munge_sumstats.py \
        --out ./Munged/$BASENAME \
        --sumstats $FILE \
        --merge-alleles eur_w_ld_chr/w_hm3.snplist \
        --N-col N \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --p Pvalue \
        --signed-sumstats Zscore,0

    # Check if munge_sumstats.py was successful
    if [ $? -ne 0 ]; then
        echo "Munge sumstats failed for $BASENAME"
        exit 1
    fi

    # Run LDSC genetic correlation
    echo "Running LDSC for $BASENAME"
    ../ldsc.py \
        --ref-ld-chr ./eur_w_ld_chr/ \
        --out ./PD_sup_results/PD_$BASENAME \
        --rg ./PD_no_UKB.sumstats.gz,./Munged/$BASENAME.sumstats.gz \
        --w-ld-chr ./eur_w_ld_chr/

    # Check if ldsc.py was successful
    if [ $? -ne 0 ]; then
        echo "LDSC failed for $BASENAME"
        exit 1
    fi
done

echo "LDSC completed"

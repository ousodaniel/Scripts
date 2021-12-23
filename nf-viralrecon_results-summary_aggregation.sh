#!/usr/bin/env bash

# Script to aggregate nf-viralrecon results (according to our setup) into one place for our downstream uses. Run from inside
# the parent run folder with the fastq dir. It takes one arg as either "ont" or "illumina".
if [[ $1 == ont ]]
then
mqc_pat=medaka
grep_pat=medaka
summ_pat=medaka

elif [[ $1 == illumina ]]
then
grep_pat=Consensus
summ_pat=variants/[bi][ov]*
mqc_pat=''
else
echo ERROR: either ont or illumina must be given as an arg; exit 1
fi

# Overwrite the dir from a previous instance of the same run of this script 
rm -r output_${PWD##*/} 2> /dev/null

mkdir  output_${PWD##*/}
agg_dir=output_${PWD##*/} && cd ${agg_dir}

mkdir -p nxt png snpEff plt dpt/amplicon dpt/genome qcs var

# Copy files of interest from various location to respective dirs
cp ../plots/amplicon/all_samples.mosdepth.heatmap.pdf plt
cp ../results/${summ_pat}/mosdepth/amplicon/*tsv dpt/amplicon
cp ../results/${summ_pat}/mosdepth/genome/*tsv dpt/genome
cp ../results/${summ_pat}/quast/transposed_report.tsv qcs
cp ../results/${summ_pat}/pangolin/*.csv png
cp ../results/${summ_pat}/snpeff/*vcf snpEff
cp ../results/multiqc/${mqc_pat}/summary_variants_metrics_mqc.csv qcs
cp ../alignment/*fasta png

# Run nexclade using the consensus from the nf-viralrecon analysis
cd nxt
sbatch -w compute06 ~/sarsis_temlate/exe/nextclade_sars.sbatch

#pid=$!

#wait 

#if [[ -e nextclade.tsv ]]; then cp nextclade.tsv nxt.tsv; fi
#cp nextclade.tsv nxt.tsv

# Merge per-sample nf-viralrecon pangolin outputs into one file
cd ../png
rm png.csv 2> /dev/null	
head_sam=`ls -1 *csv | head -n1`
head -n1 *${head_sam}* >> png-nfvr.csv

cat ./*COV*.csv | grep ${grep_pat} >> png-nfvr.csv
#if [[ `grep -v ^taxon | wc -l png-nfvr.csv` eq `ls -1 COV*.csv | wc -l` ]]; then rm COV*.csv; fi

if [[ $? == 0 ]] 
then
echo Finished result aggregation! && cd ../..
fi

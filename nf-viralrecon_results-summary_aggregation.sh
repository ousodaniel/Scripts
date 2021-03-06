#!/usr/bin/env bash

<<<<<<< HEAD
# Script to aggregate nf-viralrecon results into one place. Run from inside
# the parent run folder containing the nf-viralrecon "results" folder.
# The script is dependent on another - "nextclade_sars.sbatch" that does the
# nexclade nalysis using the consensus

if [[ $1 && $2  == ont ]]
=======
# Script to aggregate nf-viralrecon results (according to our setup) into one place for our downstream uses. Run from inside
# the parent run folder with the fastq dir. It takes one arg as either "ont" or "illumina".
if [[ $1 == ont ]]
>>>>>>> 92ff47768ed390daa0130e075002ae42d47af86f
then
mqc_pat=medaka
grep_pat=medaka
summ_pat=medaka

elif [[ $1 && $2  == illumina ]]
then
grep_pat=Consensus
summ_pat=variants/[bi][ov]* # targets the "bowtie" and "ivar" folders
mqc_pat=''
else
echo -e "ERROR: node-name and either \"ont\" or \"illumina\" must be given as args, respectively\nUSAGE: nf-viralrecon_results-summary_aggr$
fi

# Overwrite the dir from a previous instance of the same run of this script 
rm -r output_${PWD##*/} 2> /dev/null

mkdir  output_${PWD##*/}
agg_dir=output_${PWD##*/} && cd ${agg_dir}

<<<<<<< HEAD

mkdir -p nxt png snpEff plt dpt/amplicon dpt/genome qcs var

cp ../results/${summ_pat}/mosdepth/amplicon/all_samples.mosdepth.heatmap.pdf plt
=======
mkdir -p nxt png snpEff plt dpt/amplicon dpt/genome qcs var

# Copy files of interest from various location to respective dirs
cp ../plots/amplicon/all_samples.mosdepth.heatmap.pdf plt
>>>>>>> 92ff47768ed390daa0130e075002ae42d47af86f
cp ../results/${summ_pat}/mosdepth/amplicon/*tsv dpt/amplicon
cp ../results/${summ_pat}/mosdepth/genome/*tsv dpt/genome
cp ../results/${summ_pat}/quast/transposed_report.tsv qcs
cp ../results/${summ_pat}/pangolin/*.csv png
gunzip ../results/${summ_pat}/snpeff/*vcf.gz 2> /dev/null
cp ../results/${summ_pat}/snpeff/*vcf snpEff
cp ../results/multiqc/${mqc_pat}/summary_variants_metrics_mqc.csv qcs
cat ../results/${summ_pat}/*consensus* > png/${agg_dir#*_}.all.consensus.fasta

# Run nexclade using the consensus from the nf-viralrecon analysis
cd nxt
sbatch -w compute06 ~/sarsis_temlate/exe/

#pid=$!

#wait 

#if [[ -e nextclade.tsv ]]; then cp nextclade.tsv nxt.tsv; fi
#cp nextclade.tsv nxt.tsv

# Merge per-sample nf-viralrecon pangolin outputs into one file
cd ../png
rm png.csv 2> /dev/null	
head_sam=`ls -1 *csv | head -n1`
head -n1 *${head_sam}* > png-nfvr.csv

cat ./*COV*.csv | grep ${grep_pat} >> png-nfvr.csv
#if [[ `grep -v ^taxon | wc -l png-nfvr.csv` eq `ls -1 COV*.csv | wc -l` ]]; then rm COV*.csv; fi

if [[ $? == 0 ]] 
then
cd ..
echo "Finished result aggregation, see $PWD for output; the 'nxt' dir may take a bit to populate..." && cd ..
fi

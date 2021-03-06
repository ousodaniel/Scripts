#!/usr/bin/env bash

#directory hierrachy: pwd:adapter; exe; fastqc; output; primer; ref; sra
#scripts are executed from pwd

echo "##############################...Module Loading...##############################"

printf "loading modules...\n\tttrimmomatic/0.38\n\tsamtools/1.11\n\tbwa/0.7.17\n\tpython/3.7 \
\n\tbedtools/2.29.0\n\tivar/1.3\n\tsratoolkit/2.10.0\n\tnextclade/0.12.0\n"

module load trimmomatic/0.38 samtools/1.11 bwa/0.7.17 python/3.7 bedtools/2.29.0 ivar/1.3 \
nextclade/0.12.0

###################################################################################################
#retrieving slurm output file name
slurm_out=$(ls -t *.out | head -n1)

#initialise sample count
count=0

#define the files/directories for common sample in/output files of interest
#change this variables accordingly

dt=$(date '+%d-%m-%Y_%H:%M')

home_dir=~

mkdir -p ${home_dir}/variant_files_${dt} \
${home_dir}/depth_files_${dt} \
${home_dir}/nextclade_files_${dt}

var_dir=${home_dir}/variant_files_${dt}/
dph_dir=${home_dir}/depth_files_${dt}/
fastq_dir=${home_dir}/raw/batch2_01-04-2021/
nxtc_dir=${home_dir}/nextclade_files_${dt}/
touch ${dt}_nxtc.txt
nxtc_fil=${dt}_nxtc.txt
mkdir -p output/nextclade_outputs output/consensus_outputs
nxtc_out=./output/nextclade_outputs
threads=18
suff1=_1.fastq.gz
suff2=$(sed 's/1/2/' <<< ${suff1})

###################################################################################################
echo "##############################...Ref Indexing...##############################"
#### DB Indexing
#bwa: index genome
bwa index -a bwtsw -p ref/refcov ref/refcov.fa

###################################################################################################
#### Trimming(trm), Aligning(aln) & Sort(srt)
for rd1 in ${fastq_dir}*${suff1}
do
###################################################################################################
    #retrieve sample base name
    base=$(basename ${rd1} ${suff1})

    #create output directories for every sample
    mkdir -p output/${base}_outputs

    #initialise sample count
    count=$(( count + 1 ))

###################################################################################################
    #trimmomatic: trim seq adapters
    echo "##############################...Read Quality Trimming...##############################"
    trimmomatic PE -threads ${threads} ${rd1} ${fastq_dir}${base}${suff2} \
    ${base}_1.trm1.fastq ${base}_1un.trm1.fastq \
    ${base}_2.trm1.fastq ${base}_2un.trm1.fastq \
    LEADING:20 \
   MINLEN:30

###################################################################################################
    #bwmem: align reads to indexed genome (output .sam)
    #sview: convert aln-sam to aln-bam (output .bam)
    #ssort: srt aln-bam (output .bam)
    echo "##############################...Read Alignment to Ref...##############################"
    bwa mem -t ${threads} ref/refcov ${base}_1.trm1.fastq \
    ${base}_2.trm1.fastq | samtools view -hu -F 4 -F 2048 \
    | samtools sort -n -u -@ ${threads} -o ${base}_trm1.aln.srt.bam

###################################################################################################
    #sfixmate-markdup: remove duplicates
    echo "##############################...Removing Duplicates...##############################"
    time samtools fixmate -u -m ${base}_trm1.aln.srt.bam ${base}_fxm.bam
    time samtools sort -@ ${threads} -u -o ${base}_fxm.srt.bam ${base}_fxm.bam
    time samtools markdup -@ ${threads} -u -d 100 -rst -f ${base}_mkd_stats.txt ${base}_fxm.srt.bam ${base}_mkd.bam
    time samtools sort -@ ${threads} -u -o ${base}_mkd.srt.bam ${base}_mkd.bam

###################################################################################################
    #itrim: (input aln-srt.bam)trim pcr primers (requres aln-srt & primer bed)
    #(output is .bam)
    echo "##############################...iVar Primer Trim...##############################"
    ivar trim -m 20 -i ${base}_mkd.srt.bam -b primer/nCoV-2019.bed \
    -p ${base}_trm2

    #ssort: (input is aln.bam)sort trm2 bam (output is .bam)
    echo "##############################...Post iVar Sort...##############################"
    samtools sort -@ ${threads} -u -o ${base}_trm2.srt.bam ${base}_trm2.bam

###################################################################################################
    echo "##############################...Extracting Primer Read Coverage...##############################"
    #cat outputs the content of current slurm output
    #grep finds the primer coverage data (^n[C] pattern), plus its header line (-B)
    #sed replaces spaces with underscore in the header line
    cat $slurm_out | grep -B 1 '^n[C]' | sed '1 s/ /_/g' > ${base}_primers.dph

###################################################################################################

    #### ID mismatched primer to consensus
    #smpileup: (input is .bam, no ref here)pile mrg reps (output is .bam)
    #iconsensus: (input is .bam, piped) create consensus from reps mrg
    #(output is .fasta cons seq & .txt ave base Q)
    #disabled BAQ (-B) because there is no ref
    echo "##############################...Generating Consensus Genome...##############################"
    samtools mpileup -A -d 0 -Q 20 ${base}_trm2.srt.bam \
    | ivar consensus -t 0.75 -p ${base}_trm2.cnscov

    #write consensus .fasta to a text file to be used nextclade as input
    cat ${base}_trm2.cnscov.fa >> ${nxtc_fil}

    #move consensus .fasta to specific folder in the output folder
    cp ${base}_trm2.cnscov.fa output/consensus_outputs

    #bwindex: (input .fa)index consensus (output mtlp)
    echo "##############################...Indexing Reads Consensus Genome...##############################"
    bwa index -a bwtsw -p ${base}_trm2.cnscov ${base}_trm2.cnscov.fa

    #bmem: (input .fa)aln primers to consensus ref (output .sam, piped)
    #sview|ssort: convert aln-sam to aln-bam | sort aln-bam
    echo "##############################...Aligning Primers.fa to Consensus Genome...##############################"
    bwa mem -k 5 -T 16 ${base}_trm2.cnscov primer/nCoV-2019.fa \
    | samtools view -hu -F 4 | samtools sort -@ ${threads} -u -o \
    ${base}_nCoV-2019.cns.aln.srt.bam

    #smpileup: pileup primers against cns ref
    #ivariants: call primer variants from consensus (output .tsv)
    echo "##############################...Calling Primer Variants...##############################"
    samtools mpileup -A -d 0 -Q 20 --reference ${base}_trm2.cnscov.fa \
    ${base}_nCoV-2019.cns.aln.srt.bam | ivar variants \
    ${base}_nCoV-2019.cns.aln.srt.bam | ivar variants \
    -p ${base}_nCoV-2019.cns.var -t 0.25

    #get indices for mismatched primers
    echo "##############################...Converting Mismatched Primers Positions to BED...##############################"
    bedtools bamtobed -i ${base}_nCoV-2019.cns.aln.srt.bam \
    > ${base}_nCoV-2019.cns.aln.srt.bed

    echo "##############################...Getting Mismatched Primers...##############################"
    #get primers with mismatches to cns ref (output .txt)
    ivar getmasked -i ${base}_nCoV-2019.cns.var.tsv -b \
    ${base}_nCoV-2019.cns.aln.srt.bed -f \
    primer/nCoV-2019.tsv -p ${base}_nCoV-2019.msk

###################################################################################################
    echo "##############################...Trimming Mismatched-Primers Reads...##############################"
    #### Trim reads from mismatched primers
    #iremovereads: (input aln-srt-trm2(ivar trm) .bam)remove reads
    #associated with mismatched primer indeces (output .bam)*_rep${c}.trm2.srt.bam
    ivar removereads -i ${base}_trm2.srt.bam -t ${base}_nCoV-2019.msk.txt \
    -b ${base}_nCoV-2019.cns.aln.srt.bed -p ${base}_msk.trm3

    #sort trm3 bam
    echo "##############################...Sorting Alignment...##############################"
    samtools sort -@ ${threads} -u -o ${base}_msk.trm3.srt.bam ${base}_msk.trm3.bam

####################################################################################################
    echo "##############################...Calling Variants...##############################"
    #### Variant Calling
    #smpileup: pile trim3 bam
    samtools mpileup -A -d 0 -Q 20 --reference ref/refcov.fa \
    ${base}_msk.trm3.srt.bam | ivar variants -m 10 -p ${base}_variants -t 0.25 \
    -r ref/refcov.fa -g ref/refcov.gff

    #ifiltervariants: (input .tsv)filter variants (output .tsv)
    #ivar filtervariants -p ${base}_flt.variants ${base}_variants.tsv

####################################################################################################
    echo "##############################...Converting iVar .tsv to .vcf...##############################"
    ### Generating .vcf
    #ivar_variants_to_vcf.py: converting variants.tsv to variants.vcf
    python exe/ivar_variants_to_vcf.py  -po ${base}_variants.tsv ${base}_variants.vcf

####################################################################################################
    echo "##############################...Annotating Variants...##############################"
    ### Annotate variants
    #snpeff: (input .vcf) (output ann.vcf)
    snpeff="java -Xmx1g -Xmx8g -jar exe/snpEff/snpEff.jar -c exe/snpEff/snpEff.config"
    ${snpeff} covHu.1 ${base}_variants.vcf > ${base}_variants.ann.vcf

####################################################################################################
    echo "##############################...Read Depth Computation...##############################"
    #### Compute the read depths at each position for comparison
    #sdepth: (input is aln-srt.bam) depths after primer trim (output is .tsv)
    samtools depth -aa ${base}_trm2.srt.bam > ${base}_prm-trm.dph

    #sdepth: (input is aln-srt.bam) depths after masking trim (output is .tsv)
    samtools depth -aa ${base}_msk.trm3.srt.bam > ${base}_msk.dph

    #sdepth: (input is aln-srt.bam) depths after dedup trim (output is .tsv)
    samtools depth -aa ${base}_mkd.srt.bam > ${base}_mkd.dph

####################################################################################################
    echo "##############################...Archiving...##############################"
    #filter out iVar-specific warnings and store the slurm.out report for the sample
    cat $slurm_out | grep -vw '^WARNING\|^iVar\|^It' | sed -r '/^\s*$/d' > ${base}_slurm.out

    #copy variants' files to a specific folder
    cp *variants.tsv *.vcf ${var_dir}

    #copy all depth files to specific folder
    cp *.dph ${dph_dir}

    #move all output associated with every sample to a respective directory
    mv ${base}* snpEff_genes.txt snpEff_summary.html output/${base}_outputs/

    #wipe slurm.out clean
    > $slurm_out

done

echo "##############################...Clade Assignment & Tree Generation...##############################"

base=$(basename ${nxtc_fil} .txt)

nextclade.js -i ${nxtc_fil} \
--output-tsv=${base}.tsv \
--output-tree=${base}.tree.json

#copy the consolidated consensus .fasta files to a specific folder
cp ${base}* ${nxtc_dir}

#move all nextclade associated files to a designate directory
mv ${base}* ${nxtc_out}

echo "##############################...Processed ${count} Samples...##############################"



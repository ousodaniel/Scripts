set -e
echo "loading modules..."
module load trimmomatic/0.38 samtools/1.11 bwa/0.7.17 python/3.7 bedtools/2.29.0 ivar/1.3 \
sratoolkit/2.10.0

###################################################################################################
echo "dowloading reference files..."
#### Downloading files...
file="sraAccList.txt"
lines=$(cat ${file})
for l in ${lines}
do
	if [[ ! -f ${l}* ]]; then \
	echo "downloading SRA file ${l}..."; \
	fasterq-dump --split-files ${l}; fi
done

if [[ ! -f *.fa ]]; then \
echo "downloading refcov .fa..."; \
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna; \
| gunzip | mv data/refcov.fa; fi

if [[ ! -f *.fai ]]; then \
echo "downloading refcov .fai..."; \
wget https://raw.githubusercontent.com/artic-network/fieldbioinformatics/master/test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.reference.fasta.fai; \
| gunzip | mv refcov.fai; fi

if [[ ! -f *.gff ]]; then \
echo "downloading refcov .gff..."; \
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff; \
| gunzip | mv refcov.gff; fi


###################################################################################################
echo "indexing refcov genome..."
#### DB Indexing
#bwa: index genome
bwa index -a bwtsw -p refcov data/refcov.fa


###################################################################################################
echo "trimming, aligning & sorting files..."
#### Trimming(trm), Aligning(aln) & Sort(srt)
count=1
for f1 in *_1.fastq
do
	base1=$(basename ${f1} _1.fastq)
	#trimmomatic: trim seq adapters
	#trimmomatic PE ${f1} ${base1}_2.fastq \
	#${base1}_1.trm1.fastq ${base1}_1un.trm1.fastq \
	#${base1}_2.trm1.fastq ${base1}_2un.trm1.fastq \
	#ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

	#bwmem: align reads to indexed genome (output .sam)
	#sview: convert aln-sam to aln-bam (output .bam)
	#ssort: srt aln-bam (output .bam)
	bwa mem -t 4 data/db/refcov ${base1}_1.fastq \
	${base1}_2.fastq | samtools view -hb -F 4 -F 2048 \
	| samtools sort -o ${base1}_rep${count}.trm1.aln.srt.bam

	#itrim: (input aln-srt.bam)trim pcr primers (requres aln-srt & primer bed)
	#(output is .bam)
	ivar trim -i ${base1}_rep${count}.trm1.aln.srt.bam -b data/primers/nCoV-2019.bed \
	-p ${base1}_rep${count}.trm2

	
	#ssort: (input is aln.bam)sort trm2 bam (output is .bam)
	samtools sort -o ${base1}_rep${count}.trm2.srt.bam ${base1}_rep${count}.trm2.bam
	echo "replicate ${count} is sorted!"
	count=$(( count + 1 ))
done


###################################################################################################
echo "computing combined read depths before & after primer trim..."
#### Compute the read depths at each position for comparison
a= $(ls *.trm1.aln.srt.bam | wc -l)
let b="a / 2"
c=1
until [[ b -eq 0 ]]
do
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa *_rep${c}.trm1.aln.srt.bam *_rep$(( c + 1 )).trm1.aln.srt.bam > \
	Set${c}.trm1.dph
	c=$(( c + 2 ))
	c=$(( b - 1 ))
done

a= $(ls *.trm2.srt.bam | wc -l)
let b="a / 2"
c= 1
until [[ b -eq 0 ]]
do
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa *_rep${c}.trm2.srt.bam *_rep$(( c + 1 )).trm2.srt.bam > \
	Set${c}.trm2.dph
	c=$(( c + 2 ))
	c=$(( b - 1 ))
done

echo "Finished step one!\n"


###################################################################################################
echo "identifying mismatched primers..."
#### ID mismatched primer to consensus
a=$(ls *.trm2.srt.bam | wc -l)
let b="a / 2"
c=1
until [[ b -eq 0 ]]
do
	#samtools: (input is .bam)merge trm2 replicates (output is .bam)
	samtools merge Set${c}_mrg.bam *_rep${c}.trm2.srt.bam \
	*_rep$(( c + 1 )).trm2.srt.bam

	#smpileup: (input is .bam, no ref here)pile mrg reps (output is .bam)
	#iconsensus: (input is .bam, piped) create consensus from reps mrg
	#(output is .fasta cons seq & .txt ave base Q)
	samtools mpileup -A -d 0 -Q 20 Set${c}_mrg.bam \
	| ivar consensus -p Set${c}_mrg.cnscov

	#bindex: (input .fa)index consensus (output mtlp)
	bwa index -a bwtsw -p Set${c}_mrg.cnscov Set${c}_mrg.cnscov.fa
	
	#bmem: (input .fa)aln primers to consensus ref (output .sam, piped)
	#sview|ssort: convert aln-sam to aln-bam | sort aln-bam
	bwa mem -k 5 -T 16 Set${c}_mrg.cnscov data/primers/nCoV-2019.fa \
	| samtools view -hb -F 4 | samtools sort -o \
	data/primers/Set${c}_nCoV-2019.cns.aln.srt.bam
	
	#smpileup: pileup primers against cns ref
	#ivariants: call primer variants from consensus (output .tsv)
	samtools mpileup -A -d 0 -Q 20 --reference Set${c}_mrg.cnscov.fa \
	data/primers/Set${c}_nCoV-2019.cns.aln.srt.bam | ivar variants \
	-p data/primers/Set${c}_nCoV-2019.cns.var -t 0.03
	
	#get indices for mismatched primers
	bedtools bamtobed -i data/primers/Set${c}_nCoV-2019.cns.aln.srt.bam \
	> data/primers/Set${c}_nCoV-2019.cns.aln.srt.bed
	
	#get primers with mismatches to cns ref (output .txt)
	ivar getmasked -i data/primers/Set${c}_nCoV-2019.cns.var.tsv -b \
	data/primers/Set${c}_nCoV-2019.cns.aln.srt.bed -f \
	data/primers/nCoV-2019.tsv -p data/primers/Set${c}_nCoV-2019.msk
	c=$(( c + 2 ))
	b=$(( b - 1 ))
done

echo "Finished step two!\n"


###################################################################################################
echo "trimming reads with mismatched primers..."
#### Trim reads from mismatched primers
x=$(ls data/primers/*nCoV-2019.msk.txt | wc -l)
#let b="a / 2"
z=1
until [[ x -eq 0 ]]
do
	a=$(ls *.trm2.srt.bam | wc -l)
	#let b="a / 2"
	c=1
	until [[ a -eq 0 ]]
	do
		#iremovereads: (input aln-srt-trm2(ivar trm) .bam)remove reads 
		#associated with mismatched primer indeces (why not on the merged?) (output .bam)*_rep${c}.trm2.srt.bam
		f=$(echo *_rep${c}.trm2.srt.bam)
		base=$(basename ${f} _rep${c}.trm2.srt.bam)
		ivar removereads -i ${f} -t data/primers/Set${z}_nCoV-2019.msk.txt \
		-b data/primers/Set${z}_nCoV-2019.cns.aln.srt.bed \
		-p ${base}_rep${c}.msk.trm3
		#sort trm3 bam
		samtools sort -o ${base}_rep${c}.msk.trm3.srt.bam ${base}_rep${c}.msk.trm3.bam
		c=$(( c + 1 ))
		a=$(( a - 1 ))
	done
	z=$(( z + 1 ))
	x=$(( x - 1 ))
done
#get depth of trm3
a=$(ls *.msk.trm3.srt.bam | wc -l)
let b="a / 2"
c=1
until [[ b -eq 0 ]]
do
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa *_rep${c}.msk.trm3.srt.bam *_rep$(( c + 1 )).msk.trm3.srt.bam > \
	Set${c}.trm3.dph
	c=$(( c + 2 ))
	b=$(( b - 1 ))
done

echo "Finished step three!\n"

####################################################################################################
echo "calling variants..."
#### Variant Calling
a=$(ls *.msk.trm3.srt.bam | wc -l)
#let b="a / 2"
c=1
until [[ a -eq 0 ]]
do
	f7=$(echo *_rep${c}.msk.trm3.srt.bam)
	base7=$(basename ${f7} _rep${c}.msk.trm3.srt.bam)
	#smpileup: pile trim3 bam
	samtools mpileup -A -d 0 -Q 20 --reference refcov.fa \
	${f7} | ivar variants -p ${base7}_rep${c}.variants -t 0.03 \
	-r refcov.fa -g refcov.gff
	c=$(( c + 1 ))
	a=$(( a - 1 ))
done

a=$(ls *.variants.tsv | wc -l)
let b="a / 2"
c=1
until [[ b -eq 0 ]]
do
	#ifiltervariants: (input .tsv)filter variants (output .tsv)
	ivar filtervariants -p Set${c}_flt.variants *_rep${c}.variants.tsv *_rep$(( c + 1 )).variants.tsv
	c=$(( c + 2 ))
	b=$(( b - 1 ))
done

echo "Finished the pipeline!\n"
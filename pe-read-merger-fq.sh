#!/usr/bin/env bash
# for merging dual-PE (single-lane) batch data for the same sample additionally sequenced
# specify the folder containing the initial sequences (batch n-1) - "dirA"
# have the latter sequence batch (n) in a directory "dirB"

set -ex

#change the following variables accordingly
dirA=$PWD
dirB=~/raw/batch4_miseq_19-04-2021
sep=_ #this considers that the original samples name does not contain underscore(s) except such as appen$
dirOut=~/raw/batch3-4_merged
suff=_1.fastq.gz

dt=$(date +%Y-%m-%d)
log=${dt}_merge-log.txt #log file

echo -e This script does not modify your inputs whatsoever. If the output directory is identical to eith$
echo

#logging
echo MERGE_LOG: $(date '+%a %d %b %H:%M %Y') >> ${dirOut}/${log}
echo PREV_SEQ_DIR: ${dirA} >> ${dirOut}/${log}
echo CURR_SEQ_DIR: ${dirB} >> ${dirOut}/${log}
echo OUTPUT_DIR: ${dirOut} >> ${dirOut}/${log}
echo  >> ${dirOut}/${log}
echo -e OUTCOME$'\t'SAMPLE_ID$'\t'COMMENT$'\t'DATE >> ${dirOut}/${log}

count1=`ls -1 ${dirA}/*${suff} | wc -l`

if [ $count1 != 0 ] #are the files there in the first place
then
for rd in ${dirA}/*${suff} #loop through the fastq file in folder
do
  	rd1=$(basename $rd)
        rd2=$(sed 's/_1.fastq/_2.fastq/' <<< ${rd1})
        suffix=$(echo $rd1 | rev | cut -d${sep} -f1-2 | rev) # get the dynamic suffix
        base=${rd1%_[a-Z]*_[1,2].fastq.gz}
        baseHyphen=`sed 's/_/-/g' <<< ${base}`
        suffix2=`sed 's/_1.fastq.gz/_2.fastq.gz/' <<< ${suffix}` #get suffix for read 2
#	base=$(basename ${rd1} ${suffix}) # get the unique prefix
        if [[ `ls -1 ${dirB}/${base}_* 2>/dev/null | wc -l` == 2 && `ls -1 ${dirOut}/${base}_* 2>/dev/nu$
#if sample is in the latest folder of sequenced samples and not in the ouput dir
        then
            	cat ${dirA}/${rd1} ${dirB}/${base}_*1.fastq.gz > ${dirOut}/${baseHyphen}_${suffix}
                cat ${dirA}/${rd2} ${dirB}/${base}_*2.fastq.gz > ${dirOut}/${baseHyphen}_${suffix2}
                linA=$(cut -d$' ' -f1 <<< $(wc -l ${dirOut}/${baseHyphen}_${suffix})) #confirm the numbe$
                linB=$(cut -d$' ' -f1 <<< $(wc -l ${dirA}/${rd1} ${dirB}/${base}_*1.fastq.gz | tail -n1)$
                if [[ $linA == $linB ]] #test if input and output number of reads match
                then
                    	echo -e MERGED$'\t'${base}$'\t'$(date '+%a %d %b %H:%M %Y')$'\t'Reseq >> ${dirOu$
                else
                    	echo ERROR: The merge output for $base may be wrong. The number of output lines $
                        echo -e ERROR$'\t'${base}$'\t'$(date '+%a %d %b %H:%M %Y') >> ${dirOut}/${log}
                        exit 1
                fi
        elif [[ `ls -1 ${dirOut}/${base}_* 2>/dev/null | wc -l` == 2 ]]
        then
            	echo -e SKIPPED_PREV$'\t'${base}$'\t'$(date '+%a %d %b %H:%M %Y')$'\t'Already in output $
        else
            	cp ${dirA}/${base}_*1.fastq.gz ${dirOut}/${baseHyphen}_${suffix}
                cp ${dirA}/${base}_*2.fastq.gz ${dirOut}/${baseHyphen}_${suffix2}
                echo "Sample $base was not re-sequenced in the latter batch; copied it from the previous$
                echo -e COPIED-PREV$'\t'${base}$'\t'$(date '+%a %d %b %H:%M %Y')$'\t'Not reseq >> ${dirO$
        fi
done
else
    	echo "There is no file with ${suff} suffix in the directory"
fi

count2=`ls -1 ${dirB}/*${suff} 2>/dev/null | wc -l`

if [ $count2 != 0 ] # are the files remaining in the first place
then
#loop over files in     latest dir and copy if not in output dir
for rd1 in ${dirB}/*${suff}
do
  	base=$(basename ${rd1})
        prefix=${base%_[a-Z]*_[12].fastq.gz}
        prefixHyphen=`sed 's/_/-/g' <<< ${prefix}` #replace _ with -
        ex_status=`ls -1 ${dirOut}/${prefixHyphen}_* &>/dev/null; echo $?`

        if [[ $ex_status == 0 ]]
        then
            	echo -e SKIPPED_CURR$'\t'${prefix}$'\t'$(date '+%a %d %b %H:%M %Y')$'\t'Already in outpu$
        else
            	suffix=$(echo $rd1 | rev | cut -d${sep} -f1-2 | rev) #get the suffix added by the sequen$
                suffix2=`sed 's/_1.fastq.gz/_2.fastq.gz/' <<< ${suffix}` #get suffix for read 2
                cp  ${dirB}/${prefix}_*1.fastq.gz ${dirOut}/${prefixHyphen}_${suffix}
                cp  ${dirB}/${prefix}_*2.fastq.gz ${dirOut}/${prefixHyphen}_${suffix2}
                echo "Sample $prefix is newly sequenced in the latest batch; copied it from the latest b$
                echo -e COPIED-CURR$'\t'${prefix}$'\t'$(date '+%a %d %b %H:%M %Y')$'\t'New seq >> ${dirO$
        fi
done
else
    	echo #:
fi

cp ${dirOut}/${log} ${dirA}
cp ${dirOut}/${log} ${dirB}
echo -e  >> ${dirOut}/${log}
echo -e '==================================================================' >> ${dirOut}/${log}

merged=`grep -c ^MERGED ${dirOut}/${log}`
true_count=$(expr ${count1} + ${count2} - ${merged} \* 2)
if [[ `ls -1 ${dirOut}/*${suff} | wc -l` == ${true_count} ]]
then
    	echo -e The merge process finished successfully. `ls -1 ${dirOut}/*${suff} | wc -l` samples were$
else
    	echo -e Some files might have been missed in the merging process, verify your output. | tee -a $$
fi

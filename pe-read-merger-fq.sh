#!/usr/bin/env bash
# merge dual-PE (single-lane) batch data for the same sample additionally sequenced
# run script from inside the folder containing the initial sequences (batch n-1) - dirA
# have the latter sequence batch (n) in a directory labeled "dirB"

#mkdir -p merged # create a dir inside dirA to dump the merged sequences

dirA=$PWD
dirB=~/raw/batch4_19-04-2021
dirOut=~/raw/batch3-4_merged

suff=_1.fastq.gz

count=`ls -1 *${suff} | wc -l`

if [ $count != 0 ] # are the files there in the first place
then
for rd1 in *${suff}
do
	rd2=$(sed 's/_1.fastq/_2.fastq/' <<< ${rd1})
	suffix=$(echo $rd1 | rev | cut -d_ -f1-2 | rev)
  	base=$(basename ${rd1} ${suffix}) # get the unique prefix
  	if [[ `ls -1 ${dirB}/${base}* 2>/dev/null | wc -l` == 2 ]]
  	then
		cat ${rd1} ${dirB}/${base}*1* > ${dirOut}/${rd1}
		cat ${rd2} ${dirB}/${base}*2* > ${dirOut}/${rd2}
		linA=$(cut -d$' ' -f1 <<< $(wc -l ${dirOut}/${rd1}))
        	linB=$(cut -d$' ' -f1 <<< $(wc -l ${rd1} ${dirB}/${base}*1* | tail -n1))
	        if [[ $linA == $linB ]]
        	then
        		:
	    	else
    			echo ERROR: The merge output for $base may be wrong. The number of output lines do not match combined input lines 2>&1
	    		exit 1
		fi
		if [[ `ls -1 ${dirB}/${base}*1* 2>/dev/null | wc -l` == 1 ]]
		then
#			rm ${dirB}/${base}*1* ${dirB}/${base}*2*
		else
			:
	else
		echo "Sample $base was not re-sequenced in the latter batch; copied it from the previous batch's sequencing in $(basename  $dirA) to $(basename $dirOut)" > ${dirOut}/log.txt
    		cp ./${base}* ${dirOut}/
	fi
done
cp  ${dirB}/*.fastq.gz ${dirOut}/
else
    	echo "There is no file with ${suff} suffix in the directory"
fi

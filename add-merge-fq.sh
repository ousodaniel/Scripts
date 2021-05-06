#!/usr/bin/env bash
# merge dual-PE (single-lane) batch data for the same sample additionally sequenced
# run script from inside the folder containing the initial sequences (batch n-1) - dirA
# have the latter sequence batch (n) in a directory labeled "dirB" inside dirA

mkdir -p merged # create a dir inside dirA to dump the merged sequences


suff=_1.fastq.gz

count=`ls -1 *${suff} | wc -l`

if [ $count != 0 ] # are the files there in the first place
then
for rd1 in *${suff}
do
  	rd2=$(sed 's/1/2/' <<< ${rd1})
  	base=$(cut -d_ -f1 <<< ${rd1}) # get the unique cov prefix
  	count=`ls -1 dirB/${base}* | wc -l`
  	if [ $count == 2 ]
  	then
#  		rd2=$(sed 's/1/2/' <<< ${rd1})
#  		base=$(basename ${rd1} ${rmsuff})
        cat ${rd1} dirB/${base}*1* > merged/${rd1}
        cat ${rd2} dirB/${base}*2* > merged/${rd2}
        linB=$(cut -d$' ' -f1 <<< $(wc -l merged/${rd1}))
        linA=$(cut -d$' ' -f1 <<< $(wc -l ${rd1} dirB/${rd1} | tail -n1))
        if [[ $linB == $linA ]]
        then
        	:
    	else
    		echo ERROR: The merge output for $base is wrong. The number of output lines do not match combined input lines 1>&2
    		exit 1
		fi
    else
    	cp ${rd1} ${rd2} merged/
	fi
done
else
    	echo "There is no file with ${suff} suffix in the directory"
fi

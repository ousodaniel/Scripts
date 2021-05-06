#!/usr/bin/env bash
getopts ':s:' opt

if [ $opt = '-s' ]
then

	rmsuff1=$OPTARG
	rmsuff2=$(sed 's/R1/R2/' <<< ${rmsuff1})
	suff1=_1.fastq.gz
	suff2=$(sed 's/1/2/' <<< ${suff1})

	count=`ls -1 *${rmsuff1} | wc -l`

	if [ $count != 0 ]
	then
		for rd1 in *${rmsuff1}
		do
		  	base=$(basename ${rd1} ${rmsuff1})
		    mv ${rd1} ${base}${suff1}
		    mv ${base}${rmsuff2} ${base}${suff2}
		done
	else
	    echo "There is no file with ${rmsuff1} suffix in the directory"
	fi
else
	echo "Syntax error: use rename.sh -s <suffix>"
fi
#!/usr/bin/env bash
dt=$(date '+%d-%m-%Y')
mkdir -p bbduk_metrics_${dt}
met_dir=bbduk_metrics_${dt}
touch ${met_dir}/bbduk_metrics_${dt}.tsv
met_file=bbduk_metrics_${dt}.tsv

printf "file_name\tnum_reads\tnum_contam_reads\tpercnt_contam\n" > ${met_dir}/${met_file}

for file in *duk.txt
do
        file_name=$(basename ${file} .duk.txt)
        num_reads=$(grep -w '#Total' ${file} | cut -f 2)
        num_contam_reads=$(grep -w '#Matched' ${file} | cut -f 2)
        percnt_contam=$(grep -w '#Matched' ${file} | cut -f 3 | sed 's/%//g')
        echo -e "${file_name}\t${num_reads}\t${num_contam_reads}\t${percnt_contam}" >> \
        ${met_dir}/${met_file}

        #rm ${file}
done


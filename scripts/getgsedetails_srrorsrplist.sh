#!/bin/bash
# Author: A.T
# Need efetch/esearch from entrez-direct. Install via conda

#set -euo pipefail
if [ "$#" -ne 3 ]; then
echo "bash script.sh srr/srp file(one srr/srp per line) SRR/SRP outputname"
echo "bash script.sh srrid.txt SRR myoutput"
exit
fi

# check for input files
if [ ! -f $1 ];then
        echo "Inputfile:${1} not found! Exiting"
        exit
fi


id=$(echo "$2" | awk '{print tolower($0)}')
#echo "id:${id}"
if [[ $id == "srr" ]];then
        grep "^SRR" ${1} > ${1}.temp
elif [[ $id == "srp" ]]; then
        grep "^SRP" ${1} > ${1}.temp
else
        echo "Invalid argument 2. Only SRR or SRP"
        exit
fi
path=$PWD

total=$(wc -l ${1}|cut -d' ' -f1)
srrtotal=$(wc -l ${1}.temp|cut -d' ' -f1)
echo "Valid # of SRR/SPRs=${srrtotal} out of ${total} ids given"
tempdir=$(mktemp -d)
echo "tempdir:${tempdir}"
cd $tempdir
split -l 1  ${path}/${1}.temp
echo "Finding SRR for each SRPs"
ls -1 x*|while read i;do
        cat $i|while read j;do
                srr=$j
                esearch -db sra -query $srr|efetch -format runinfo|sed '1d'|grep $srr |cut -d',' -f 1,30,21|sed '/^$/d'|tr "," "\t"
                sleep 3
                done
        done >  ${path}/${3}.srr_srp_gsmid.txt
rm -rf ${path}/${1}.temp
echo "SRR SRP GSM list created"
rscriptpath=get_srr-srp-gse_details_tab.R
if [ ! -f ${rscriptpath} ]; then
        echo "Rscript ${rscriptpath} not found!"
        echo "Needs ${rscriptpath} in path! Exiting!"
        exit
fi
cd ${path}
Rscript ${rscriptpath} ${3}.srr_srp_gsmid.txt $3

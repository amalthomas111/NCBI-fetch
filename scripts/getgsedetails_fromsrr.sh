#!/bin/bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
echo "sh script.sh srrfile(one srr per line) outputname"
exit
fi
path=$PWD
grep "^SRR" ${1} > ${1}.temp
total=$(wc -l ${1}|cut -d' ' -f1)
srrtotal=$(wc -l ${1}.temp|cut -d' ' -f1)
echo "# of SRRs=${srrtotal} out of ${total} ids given"
tempdir=$(mktemp -d)
echo "tempdir:${tempdir}"
cd $tempdir
split -l 1  ${path}/${1}.temp
ls -1 x*|while read i;do cat $i|while read j;do srr=$j; esearch -db sra -query $srr|efetch -format runinfo|sed '1d'|grep $srr |cut -d',' -f 1,30,21|sed '/^$/d'|tr "," "\t"  >> ${path}/${2}.srr_srp_gsmid.txt;done;done
cd ${path}
Rscript ~/panases_soft/myexec/get_srr-srp-gse_details_tab.R ${2}.srr_srp_gsmid.txt $2
rm -rf ${1}.temp

#!/bin/bash
# Author: A.T
# Needs srapath from NCBI toolkit, active internet connection
# This script finds  new 10x SRPs and SRR that are not processed.
# Input file is the ncbi_fetch outputfile for 10x
# Find SRP with 10x bam file in GEO. For those SRPs without bam files, write
# to separate file so that it can be processed from fastqs

#set -euo pipefail
if [ "$#" -ne 1 ]; then
echo "bash script.sh <ncbi_fetch10xoutputfile>"
exit
fi

#ncbi_fetch_outputfile
#ncbi_fetch="Geo_10x_dataset.tsv"
ncbi_fetch=$1

#finished human srp file
human_srp="finished_human_srps.txt"
human_analysis_path_HPC=/home/cmb-06/as/desenabr/scream/db/hg38/analyses

if [ ! -f ${human_srp} ]; then
        echo "Finished human srp found! Creating SRP from HPC path"
        ls -1 ${human_analysis_path_HPC}|cut -d_ -f1 | sort -u | grep SRP > \
                finished_human_srps.txt
fi

#finished mouse srp file
mouse_srp="finished_mouse_srps.txt"
mouse_analysis_path_HPC=/home/cmb-06/as/desenabr/scream/db/mm10/analyses
if [ ! -f ${mouse_srp} ]; then
        echo "Finished mouse srp found! Creating SRP from HPC path"
        ls -1 ${mouse_analysis_path_HPC}|cut -d_ -f1 | sort -u | grep SRP > \
                finished_mouse_srps.txt
fi

#finished humanmouse srp file
humanmouse_srp="finished_humanmouse_srps.txt"
cat ${human_srp} ${mouse_srp} |sort -u > finished_humanmouse_srps.txt


date=$(date +%Y-%m-%d)
tempdir="other"
totaldir="total"

if [ -f ${ncbi_fetch} ]; then
mkdir -p ${tempdir}
mkdir -p ${totaldir}
grep "Mus musculus" ${ncbi_fetch} |grep -v "Homo sapiens" > \
                               ${tempdir}/Geo_10x_dataset_onlyMouse.tsv
grep "Homo sapiens" ${ncbi_fetch} |grep -v "Mus musculus" > \
                               ${tempdir}/Geo_10x_dataset_onlyHuman.tsv
grep "Homo sapiens" ${ncbi_fetch} |grep "Mus musculus" > \
                                ${tempdir}/Geo_10x_dataset_humanmouse.tsv
else
        echo "ncbi fetch outputfile (Geo_10x_dataset.tsv) not found!"
        echo "First run python3 getGSE_fromNCBI_10x.py. Exiting!" 
        exit
fi

#get the new srps
cut -f3 ${tempdir}/Geo_10x_dataset_onlyHuman.tsv|sed '/^$/d' |sort -u > \
                                ${totaldir}/total_human_SRPs_${date}.txt
cut -f3 ${tempdir}/Geo_10x_dataset_onlyMouse.tsv|sed '/^$/d' |sort -u > \
                                ${totaldir}/total_mouse_SRPs_${date}.txt
cut -f3 ${tempdir}/Geo_10x_dataset_humanmouse.tsv|sed '/^$/d' |sort -u > \
                                ${totaldir}/total_humanmouse_SRPs_${date}.txt

if [ -f ${human_srp} ]; then
        echo "${human_srp} file found!"
        grep -v -f ${human_srp} ${totaldir}/total_human_SRPs_${date}.txt > \
                                ${totaldir}/new_human_10xsrps_${date}.txt
        grep -f ${totaldir}/new_human_10xsrps_${date}.txt ${tempdir}/Geo_10x_dataset_onlyHuman.tsv|\
                awk  'BEGIN{FS="\t"}{if($6<=10)print $0}' > \
                ${tempdir}/new_10xhuman_details_max10GSMs_${date}.tsv
        cut -f3 ${tempdir}/new_10xhuman_details_max10GSMs_${date}.tsv |sed '/^$/d'|\
                sort -u > ${tempdir}/new_10xhuman_SRPs_withmax10GSMsperSRP_${date}.txt
        cat ${tempdir}/new_10xhuman_SRPs_withmax10GSMsperSRP_${date}.txt |\
                while read i;do
                        url=$(srapath $i -f names --raw -p typ=srapub_files  | \
                        grep srapub_files | cut '-d|' -f8|sed '/^$/d'| \
                        grep -e "bam" -e "BAM"|grep -v "Cochlea")
                        if [ -z "$url" ]; then
                                echo -e ${i} >> ${tempdir}/newhuman_10xbam_SRPsnobam.txt
                        else
                                for j in $url;do
                                        echo -e  $i"\t"$j >> ${tempdir}/newhuman_10xbam_url.txt
                                done
                        fi
                        done
        cat ${tempdir}/newhuman_10xbam_SRPsnobam.txt|sort -u > newhuman_10xbam_SRPsnobam_${date}.txt
        cat ${tempdir}/newhuman_10xbam_url.txt|sort -u > newhuman_10xbam_url_${date}.txt
        echo "human done"
else
        echo "${human_srp} file not found"
fi

if [ -f ${mouse_srp} ]; then
        echo "${mouse_srp} file found!"
        grep -v -f ${mouse_srp} ${totaldir}/total_mouse_SRPs_${date}.txt > \
                ${totaldir}/new_mouse_10xsrps_${date}.txt
        grep -f ${totaldir}/new_mouse_10xsrps_${date}.txt  ${tempdir}/Geo_10x_dataset_onlyMouse.tsv |\
                awk  'BEGIN{FS="\t"}{if($6<=10)print $0}' > \
                ${tempdir}/new_10xmouse_details_max10GSMs_${date}.tsv
        cut -f3 ${tempdir}/new_10xmouse_details_max10GSMs_${date}.tsv |\
                sort -u > ${tempdir}/new_10xmouse_SRPs_withmax10GSMsperSRP_${date}.txt
        cat ${tempdir}/new_10xmouse_SRPs_withmax10GSMsperSRP_${date}.txt|\
                while read i;do
                        url=$(srapath $i -f names --raw -p typ=srapub_files  | \
                        grep srapub_files | cut '-d|' -f8|sed '/^$/d'| \
                        grep -e "bam" -e "BAM"|grep -v "Cochlea")
                        if [ -z "$url" ]; then
                                echo -e ${i} >> ${tempdir}/newmouse_10xbam_SRPsnobam.txt
                        else
                                for j in $url;do
                                        echo -e  $i"\t"$j >> ${tempdir}/newmouse_10xbam_url.txt
                                done
                        fi
                        done
        cat ${tempdir}/newmouse_10xbam_SRPsnobam.txt|sort -u > newmouse_10xbam_SRPsnobam_${date}.txt
        cat ${tempdir}/newmouse_10xbam_url.txt|sort -u > newmouse_10xbam_url_${date}.txt
        echo "mouse done"
else
        echo "${mouse_srp} file not found!"
fi

if [ -f ${humanmouse_srp} ]; then
        echo "${humanmouse_srp} file found!"
        grep -v -f ${humanmouse_srp} ${totaldir}/total_humanmouse_SRPs_${date}.txt >\
                ${totaldir}/new_humanmouse_10xsrps_${date}.txt
        grep -f ${totaldir}/new_humanmouse_10xsrps_${date}.txt  ${tempdir}/Geo_10x_dataset_humanmouse.tsv|\
                awk  'BEGIN{FS="\t"}{if($6<=10)print $0}'> \
                ${tempdir}/new_10xhumanmouse_details_max10GSMs_${date}.tsv
        cut -f3 ${tempdir}/new_10xhumanmouse_details_max10GSMs_${date}.tsv|\
                sort -u > ${tempdir}/new_10xhumanmouse_SRPs_withmax10GSMsperSRP_${date}.txt
        cat ${tempdir}/new_10xhumanmouse_SRPs_withmax10GSMsperSRP_${date}.txt|\
                while read i;do
                        url=$(srapath $i -f names --raw -p typ=srapub_files  | \
                        grep srapub_files | cut '-d|' -f8|sed '/^$/d'| \
                        grep -e "bam" -e "BAM"|grep -v "Cochlea")
                        if [ -z "$url" ]; then
                                echo -e ${i} >> ${tempdir}/newhumanmouse_10xbam_SRPsnobam.txt
                        else
                                for j in $url;do
                                        echo -e  $i"\t"$j >> ${tempdir}/newhumanmouse_10xbam_url.txt
                                done
                        fi
                        done
        cat ${tempdir}/newhumanmouse_10xbam_SRPsnobam.txt|sort -u > newhumanmouse_10xbam_SRPsnobam_${date}.txt
        cat ${tempdir}/newhumanmouse_10xbam_url.txt|sort -u > newhumanmouse_10xbam_url_${date}.txt
        echo "human mouse done"
else
        echo "${humanmouse_srp} file not found!"
fi

####### getting srr/gsm for new human SRPs with no bam ###############

human_output=newhuman_10xbam_SRPsnobam_${date}.txt
mouse_output=newmouse_10xbam_SRPsnobam_${date}.txt
humanmouse_output=newhumanmouse_10xbam_SRPsnobam_${date}.txt
script=/home/rcf-40/amalthom/panfs/software/myexec/getgsedetails_srrorsrplist.sh
if [ ! -f ${script} ]; then
        echo "bash script for SRR/GSM details not found!Exiting!"
        exit
fi

if [ -f ${human_output} ]; then
        echo "Getting SRR/GSM details for human SRPs with no bam"
        bash ${script} ${human_output} SRP ${human_output/.txt/}
fi

if [ -f ${mouse_output} ]; then
        echo "Getting SRR/GSM details for mouse SRPs with no bam"
        bash ${script} ${mouse_output} SRP ${mouse_output/.txt/}
fi
if [ -f ${humanmouse_output} ]; then
        echo "Getting SRR/GSM details for humanmouse SRPs with no bam"
        bash ${script} ${humanmouse_output} SRP ${humanmouse_output/.txt/}
fi

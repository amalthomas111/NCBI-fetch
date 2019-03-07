grep "Mus musculus" Geo_10x_dataset.tsv |grep -v "Homo sapiens" > Geo_10x_dataset_onlyMouse.tsv
grep "Homo sapiens" Geo_10x_dataset.tsv |grep -v "Mus musculus" > Geo_10x_dataset_onlyHuman.tsv
grep "Homo sapiens" Geo_10x_dataset.tsv |grep "Mus musculus" > Geo_10x_dataset_humanmouse.tsv
cut -f3 Geo_10x_dataset_onlyHuman.tsv|sed '/^$/d' |sort -u > total_human_SRPs.txt
cut -f3 Geo_10x_dataset_onlyMouse.tsv|sed '/^$/d' |sort -u > total_mouse_SRPs.txt
grep -v -f finished_human_srps.txt total_human_SRPs.txt > new_human_srps.txt
grep -v -f finished_mouse_srps.txt total_mouse_SRPs.txt > new_mouse_srps.txt
mkdir -p selected;grep -f new_human_srps.txt Geo_10x_dataset_onlyHuman.tsv|awk  'BEGIN{FS="\t"}{if($6<=10)print $0}' > selected/new_10x_human_details.tsv
cut -f3 selected/new_10x_human_details.tsv |sort -u > selected/unfinished_$(date +%Y-%m-%d).onlyhuman_onlywithSRPs_max10GSMsperSRP.tsv
mkdir -p selected;grep -f new_mouse_srps.txt  Geo_10x_dataset_onlyMouse.tsv |awk  'BEGIN{FS="\t"}{if($6<=10)print $0}' > selected/new_10x_mouse_details.tsv
cut -f3 selected/new_10x_mouse_details.tsv |sort -u > selected/unfinished_$(date +%Y-%m-%d).onlymouse_onlywithSRPs_max10GSMsperSRP.tsv
cd selected;
cat unfinished_2019-03-06.onlyhuman_onlywithSRPs_max10GSMsperSRP.tsv |while read i;do srapath $i -f names --raw -p typ=srapub_files  | grep srapub_files | cut '-d|' -f8|sed '/^$/d'|grep -e "bam" -e "BAM"|grep -v "Cochlea";done > newhuman_bam_url.txt
cat unfinished_2019-03-06.onlymouse_onlywithSRPs_max10GSMsperSRP.tsv|while read i;do srapath $i -f names --raw -p typ=srapub_files  | grep srapub_files | cut '-d|' -f8|sed '/^$/d'|grep -e "bam" -e "BAM"|grep -v "Cochlea";done > newmouse_bam_url.txt

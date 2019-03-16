# NCBI-fetch
Scripts to query &amp; search data from GEO

getGSE_fromNCBI_10x.py: python3 script to search in NCBI GEO for keywords related to single-cell RNAseq.
Dependencies: helper.py
Generates the following outputs
 - Geo_sc_Datasets.tsv : Information about studies that are found
 - log.txt: stdout &amp stderr streams
 - sc_protocol_complete.tsv: Information about library preparation & growth protocols
 - Geo_10x_dataset.tsv : studies that from 10x chromium platform
 - uid_finished.txt: Query uids that are arleady fetched

getgsedetails_fromsrr.sh: BASH script to get the metadata related to a SRR id. Needs get_srr-srp-gse_details_tab.R script.
get10xnewLibs.sh takes the output files and find new 10x studies that are not present in our collection

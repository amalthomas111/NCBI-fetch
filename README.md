# NCBI-fetch
Scripts to query &amp; search data from GEO

getGSE_fromNCBI_10x.py searches GEO for keywords related to single-cell RNAseq and outputs
 - Geo_sc_Datasets.tsv : Information about studies that are found
 - sc_protocol_complete.tsv: Information about library preparation & growth protocols
 - Geo_10x_dataset.tsv : studies that from 10x chromium platform
 - uid_finished.txt: Query uids that are arleady fetched
 - get10xnewLibs.sh takes the output files and find new 10x studies that are not present in our collection

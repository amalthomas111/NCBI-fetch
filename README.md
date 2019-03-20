# NCBI-fetch
Scripts to query &amp; search data from GEO

getGSE_fromNCBI_10x.py: python3 script to search in NCBI GEO for keywords 
mentioned in the input JSON file.

Requirements:
- Needs **active internet connection**
- JSON input file: The script requires an input JSON file with query keywords (s) mentioned as a list.
- `helper.py`: Some helper class/functions for the script are declared in `helper.py` script.

To run:
```
python3 getGSE_fromNCBI_10x.py <input json file>
```


Generates the following outputs
 - Geo_sc_Datasets.tsv : Information about studies that are found
 - log.txt: STDOUT &amp; STDERR streams
 - sc_protocol_complete.tsv: Information about library preparation & growth protocols
 - Geo_10x_dataset.tsv : studies that from 10x chromium platform
 - uid_finished.txt: Query uids that are arleady fetched

`getgsedetails_srrorsrplist.sh`: BASH script to get the metadata related to a SRR id or SRP id. 

Requirements:
- Needs **active internet connection**
- Needs efetch/esearch from entrez-direct. Install via conda (```conda install -c bioconda entrez-direct```)
- Needs get_srr-srp-gse_details_tab.R script in the path. This Rscript needs GEOquery &amp; data.table packages
- SRR id or SRP id input file (one per line).

An e.g. command to get metadata for SRR ids is:
```
bash getgsedetails_srrorsrplist.sh srrid.txt SRR mysrrinformation
```

An e.g. command to get metadata for SRP ids is:
```
bash getgsedetails_srrorsrplist.sh srpids.txt SRP mysrrinformation
```
<!---
get10xnewLibs.sh takes the output files and find new 10x studies that are not present in our collection
-->

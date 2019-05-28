#' ---
#' Rscript to "Query & get # of hits in NCBI db such as pubmed/gds"
#' author: A.T
#' To run: RScript get_noofhits_NCBI.R
#' ---
#' 
## ------------------------------------------------------------------------
suppressPackageStartupMessages({
library(RISmed)
library(ggplot2) 
})


#' 
#' 
#' # function to search db
## ------------------------------------------------------------------------
findhits = function(keyword,db="pubmed",outputname,dttype="edat"){
tally = array()
r = list()
x = 1
for (i in 2008:2019){
  Sys.sleep(1)
  ids=c()
  for(k in 1:length(keyword)){
    r[[k]] = EUtilsSummary(keyword[k], type='esearch', db=db, 
                           mindate=paste0(i,"/1/1"), maxdate=paste0(i,"/12/31"),
                           datetype=dttype)
    ids = c(ids,QueryId(r[[k]]))
  }
  tally[x] = length(unique(ids))
  x = x + 1
}
names(tally) = 2008:2019

df = data.frame(year = names(tally),f=tally)
ggplot(data=df, aes(x=year,y=f)) +
geom_bar(stat = "identity") +
geom_text(aes(label=f), vjust=-0.5) + theme_bw() +
  theme(panel.border= element_rect(colour = "black",size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")+ ylab("") + ggtitle(paste0("# of hits in ",db," with ",outputname))
ggsave(filename = paste0(db,"_",outputname,".pdf"),dpi = 300)

return(tally)
}

#' # indrop
## ------------------------------------------------------------------------
cat("\nIndrop hits:\n")
findhits(c("indrop","inDrop"),db="gds","indrop","pdat")

#' # 10x
## ------------------------------------------------------------------------
cat("\n10x hits:\n")
findhits(c('"10x chromium"','"10x Genomics"'),db="gds","10x","pdat")


#' # dropseq but not 10x no indrop
## ------------------------------------------------------------------------
outputname = "dropseq_not10xindrop"
tally = array()
r = list()
x = 1
for (i in 2008:2019){
  Sys.sleep(1)
  res1 = EUtilsSummary('dropseq', type='esearch', db='gds',
                        mindate=paste0(i,"/1/1"),
                        maxdate=paste0(i,"/12/31"),datetype="pdat")
  res2 = EUtilsSummary('drop-seq', type='esearch', db='gds',
                        mindate=paste0(i,"/1/1"),
                        maxdate=paste0(i,"/12/31"),datetype="pdat")
  totaldropseqids = unique(c(QueryId(res1),QueryId(res2)))
    
  res3 = EUtilsSummary('"10x chromium"', type='esearch', 
                        db='gds',mindate=paste0(i,"/1/1"),
                        maxdate=paste0(i,"/12/31"),datetype="pdat")
  res4 = EUtilsSummary('"10x Genomics"', type='esearch', db='gds',
                        mindate=paste0(i,"/1/1"),
                        maxdate=paste0(i,"/12/31"),datetype="pdat")
  res5 = EUtilsSummary('indrop', type='esearch', db='gds',
                       mindate=paste0(i,"/1/1"),maxdate=paste0(i,"/12/31"),
                       datetype="pdat")
  #total10x = unique(c(QueryId(res3),QueryId(res4)))
  totalother = unique(c(QueryId(res3),QueryId(res4),QueryId(res5)))
  
  uniquedropseq = totaldropseqids[!totaldropseqids %in% totalother]
  tally[x] = length(unique(uniquedropseq))
  x = x + 1
}

names(tally) = 2008:2019

df = data.frame(year = names(tally),f=tally)
ggplot(data=df, aes(x=year,y=f)) +
geom_bar(stat = "identity") +
geom_text(aes(label=f), vjust=-0.5) + theme_bw() +
  theme(panel.border= element_rect(colour = "black",size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")+ ylab("") + ggtitle(paste0("# of hits in gds with ",outputname))
ggsave(filename = paste0("gds_",outputname,".pdf"),dpi = 300)

cat("\nDropseq hits that are not 10x and indrop:\n")
tally

#' # singlecell
#' 
#' For single cell it is better to search in pubmed rather than gds
## ------------------------------------------------------------------------
cat("\nSingle-cell hits in gds:\n")
findhits(c('"single-cell RNA-Seq"','"single-cell RNASeq"',
           '"singlecell RNASeq"'),db="gds","singlecell_gds","pdat")

#' single cell search in pubmed
## ------------------------------------------------------------------------
cat("\nSingle-cell hits in pubmed:\n")
findhits(c('"single-cell RNA-Seq"','"single-cell RNASeq"',
           '"singlecell RNASeq"'),db="pubmed","singlecell_pubmed")

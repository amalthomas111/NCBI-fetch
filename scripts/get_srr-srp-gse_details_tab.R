args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript script.R <srr-srp-gse id file>  <outputname>", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input file not found!\
Usage: Rscript script.R <srr-srp-gse id file> <outputname>", call.=FALSE)
}


suppressPackageStartupMessages({
  require(GEOquery)
  require(data.table)
})

getGSEname = function(srr,gse,srp){
  tryCatch({
      gds <- getGEO(gse)
      df = data.frame(srr=srr,
                      gsm=gse,
                      srp=srp,
                      gse=Meta(gds)$series_id,
                      title=Meta(gds)$title,
                      source_name=Meta(gds)$source_name_ch1,
                      #gender= gsub(' ',"",tstrsplit(Meta(gds)$characteristics_ch1[[1]],":")[[2]]),
                      charecteristics=paste(Meta(gds)$characteristics_ch1,collapse = ","),
                      #description = Meta(gds)$description,
                      species = Meta(gds)$organism_ch1
                      )
      return(df)
  },
  error = function(e){cat("Invalid GSE:",gse,"\n")})
}

#l = c("GSM2179767","GSM2258360")
l = read.table(args[1],col.names = c("srr","srp","gsm"),stringsAsFactors = F,sep="\t")
head(l)

d = data.frame(srr=NULL,gsm=NULL,srp=NULL,gse=NULL,title=NULL,source_name=NULL,charecteristics=NULL,
               #description=NULL,
               species=NULL)
for(i in 1:nrow(l)){
  srr=l[i,1]
  gse=l[i,3]
  srp=l[i,2]
  cat(srr,gse,srp,"\n")
  d = rbind(d,getGSEname(srr,gse,srp))
}
d = d[!duplicated(d$srr),]
write.table(d,file = paste0(args[2],"_srrgsesrp_details.tsv"),sep = "\t",quote = F,row.names = F)

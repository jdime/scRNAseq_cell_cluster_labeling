####################################
### Script 'obtains_GSVA_for_MatrixColumns.R' obtains Gene Set Variation Analysis enrichment scores
### for each column of --infile_mat vs. each gene set from --infile_gmt
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_GSVA_for_MatrixColumns.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
### 'optparse'   (CRAN) to handle one-line-commands
### 'parallel'   (CRAN) to run parallel processes
### 'data.table' (CRAN) to speed up reading matrices with 'fread' instead of 'read.table'
### 'GSA'        (CRAN) to handle *gmt infile
### 'GSVA'       (bioconductor) to run the gsva function
### 'qvalue'     (bioconductor) to get FDR/q-values from GSVA's p-values
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSA))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(qvalue))
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

####################################
### Get inputs from command line argumets
####################################

option_list <- list(
  make_option(c("-i", "--infile_mat"), default="NA",
              help="A path/name to a <tab> delimited *file* with genes in rows and arrays (e.g. clusters or conditions) in columns, like:
              genes  clust1  clust2  clust3 ... etc
              RP113  0       0.0045  0.0008
              FAM14  0.0077  0.0175  0.0082
              NOC2L  0.0800  0.1532  0.0745
              ...etc"),

  make_option(c("-c", "--infile_gmt"), default="NA",
              help="A path/name to a <tab> delimited *file* of gene sets in *gmt format, like:
              GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 Gene3
              GeneSet2_ID  GeneSet2_Name  Gene4 Gene5
              ... etc"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID"),
  
  make_option(c("-e", "--pvalue_cutoff"), default="0.05",
              help="This script produce a *filtered.tsv matrix with pairs passing -e and -f filters and unfiltered outfiles"),
  
  make_option(c("-f", "--fdr_cutoff"), default="0.1",
              help="Same as -e option, but for FDR scores")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileGmt       <- opt$infile_gmt
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
PvalueCutoff    <- as.numeric(opt$pvalue_cutoff)
FdrCutoff       <- as.numeric(opt$fdr_cutoff)
Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_mat", "infile_gmt", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs and define outfiles
####################################
writeLines("\n*** Create outdirs ***\n")
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Outdir)
#
dir.create(file.path(Outdir, "GSVA"), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),        showWarnings = F, recursive = T)

OutfileEnrichmentScores<-paste(Tempdir,"/",PrefixOutfiles,".GSVA_enrichment_scores.tsv", sep="")
OutfilePvalues         <-paste(Tempdir,"/",PrefixOutfiles,".GSVA_pvalues.tsv", sep="")
OutfileFdrvalues       <-paste(Tempdir,"/",PrefixOutfiles,".GSVA_fdr_values.tsv", sep="")
OutfileAllScores       <-paste(Tempdir,"/",PrefixOutfiles,".GSVA_all_scores_table.tsv", sep="")
OutfileFilteredES      <-paste(Tempdir,"/",PrefixOutfiles,".GSVA_filtered.tsv", sep="")

####################################
### Load data
####################################
writeLines("\n*** Load data ***\n")

### Creating object fullmat with the matrix of genes(rows) vs. arrays(columns)
fullmat<-as.matrix(data.frame(fread(InfileMat, sep="\t", na.strings=c("NA")), row.names=1))

### Creates object gmt2 with the gene set memberships
gmt1<-GSA.read.gmt(InfileGmt)
gmt2<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

### Get and print out GSVA enrichment scores
StartTimeGsva<-Sys.time()
EnrichmentScores<-gsva(expr=fullmat, gset.idx.list=gmt2, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=T, parallel.sz=0)
SortedRowNames<-rownames(EnrichmentScores)
SortedColNames<-colnames(EnrichmentScores)
#
EnrichmentScores<-EnrichmentScores[SortedRowNames,SortedColNames]
write.table(data.frame("ENRICHMENT"=rownames(EnrichmentScores), EnrichmentScores), OutfileEnrichmentScores, row.names = F,sep="\t",quote = F)

### maps gene set members from gmt2
Classes.list<-NULL
for (i in 1:nrow(EnrichmentScores)){
  tmp1<-unlist(strsplit(gmt1[[2]][match(rownames(EnrichmentScores)[i],gmt1[[3]])],"%"))[3]
  tmp2<-paste(unlist(gmt2[[match(rownames(EnrichmentScores)[i],names(gmt2))]]),collapse = ",")
  Classes.list<-rbind(Classes.list,c(tmp1,tmp2))
}
EndTimeGsva<-Sys.time()

####################################
### Get and print Enrichment, P-values and Q-values (FDR)
####################################
writeLines("\n*** Get and print Enrichment, P-values and Q-values (FDR) ***\n")

### Originally used something like:
### qvalues<-qvalue(pvalues,lambda=seq(0.05,0.45,0.01))$lfdr
### But this was producing errors for some clusters given their p-values distributions, like:
### "Error in smooth.spline(lambda, pi0, df = smooth.df) : 
### missing or infinite values in inputs are not allowed"
### See https://github.com/StoreyLab/qvalue/issues/11 and 
### http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
### Javier Diaz replaced this by:
### qvalues<-qvalue(pvalues,pi0=1)$lfdr

HeaderCutoff<-paste(c("PassCutoff", "_p", PvalueCutoff, "_fdr", FdrCutoff) , sep="",collapse = "")
HeaderCutoff
HeadersForPandQvalues<-paste("CLASS","ColumnHeader","EnrichmentScore","p.Val","FDR", HeaderCutoff ,sep="\t",collapse = "")
write(x=HeadersForPandQvalues,file=OutfileAllScores)

for (columnNumber in 1:ncol(EnrichmentScores)){
  pvalues<-pnorm(-abs(scale(EnrichmentScores[,columnNumber])[,1]))
  qvalues<-qvalue(pvalues,pi0=1)$lfdr
  PassCutoffs<-ifelse((pvalues<=PvalueCutoff & qvalues<=FdrCutoff)==TRUE,1,0)
  concatenatedResults<-cbind(colnames(fullmat)[columnNumber], EnrichmentScores[,columnNumber], pvalues,qvalues,PassCutoffs)
  # Write out Table with CLASS ColumnHeader EnrichmentScore p.Val FDR PassCutoff
  write.table(x=data.frame("CLASS"=rownames(concatenatedResults),concatenatedResults),file=OutfileAllScores, row.names = F,sep="\t",quote = F,col.names = F,append = T)
}

### Index and print out GSVA p.Val and FDR
dataEPQ <- read.table(file = OutfileAllScores,row.names = NULL ,header = T)
#
ForPvaluesMat<- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"p.Val"])
PvaluesMat<-xtabs(z~x+y, data=ForPvaluesMat)
PvaluesMat<-PvaluesMat[SortedRowNames,SortedColNames]
HeadersForPvalues<-paste(c("PVALUES", colnames(PvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForPvalues,file=OutfilePvalues)
write.table(x=PvaluesMat, file=OutfilePvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
ForFDRvaluesMat <- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"FDR"])
FdrvaluesMat<-xtabs(z~x+y, data=ForFDRvaluesMat)
FdrvaluesMat<-FdrvaluesMat[SortedRowNames,SortedColNames]
HeadersForFdrvalues<-paste(c("FDR_VALUES", colnames(FdrvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForFdrvalues,file=OutfileFdrvalues)
write.table(x=FdrvaluesMat, file=OutfileFdrvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
#
FilteredESMatLogical<-(PvaluesMat<=PvalueCutoff & FdrvaluesMat<=FdrCutoff)
FilteredESMatLogical<-FilteredESMatLogical[SortedRowNames,SortedColNames]
FilteredESMatValues<-ifelse(FilteredESMatLogical==TRUE,EnrichmentScores,NA)
write.table(data.frame("ENRICHMENT_FILTERED"=rownames(EnrichmentScores),FilteredESMatValues),OutfileFilteredES, row.names = F,sep="\t",quote = F)

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".GSVA_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
writeLines("\n*** Report time used ***\n")

EndTimeOverall<-Sys.time()

TookTimeGsva    <-format(difftime(EndTimeGsva,    StartTimeGsva,    units = "secs"))
TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "secs"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".GSVA_CPUusage.txt", sep="")
ReportTime<-c(
  paste("gsva",   TookTimeGsva,   collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/GSVA/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".GSVA_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/GSVA/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

####################################
### Turning warnings on
####################################
options(warn = oldw)

####################################
### Finish
####################################

print("END - All done!!! Took time:")
print(ReportTime)

quit()

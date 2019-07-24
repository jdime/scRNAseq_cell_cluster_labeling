####################################
### Script 'obtains_METANEIGHBOR_for_MatrixColumns.R' MetaNeighborUS (unsupervised) scores quantifying cell type replicability across datasets
### using neighbor voting for each column of --infile_mat vs. each gene set from --infile_gmt
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_METANEIGHBOR_for_MatrixColumns.R -h'
### for help
####################################

####################################
### MODIFICATION TO METANEIGHBOR SOURCE CODE TO AVOID AVERAGING NEIGHBOR AUROC's
### MetaNeighborUS gets the average AUROC for each pair of cell types. It creates a symetric matrix including 'training' and 'testing' cell types
### (called studies 1 and 2) in the MetaNeighbor R vignette. For details see:
### https://bioconductor.org/packages/release/bioc/vignettes/MetaNeighbor/inst/doc/MetaNeighbor.pdf
###
### For the purpose of evaluating cell cluster labeling methods, instead of studies 1 and 2, Diaz-Mejia et al (2019)
### used cell clusters from scRNA-seq as 'testing' (i.e. study 1) cell types and cell type signatures as 'training' (i.e. study 2) cell types.
###
### The typical output of the function `MetaNeighborUS` in the MetaNeighbor source code `~/.../R/MetaNeighborUS.R`
### is a cell type-by-cell type mean AUROC matrix, which is built by treating each pair
### of cell types as both 'testing' and 'training' data for MetaNeighbor, then taking the *average* AUROC for each pair of
### neighbor AUROC scores across 'testing' and 'training'.
### The 'testing' and 'training' folds will not be identical because each test cell type is scored out of its own dataset,
### and differences in dataset heterogeneity influence scores).
### See Metaneighbor's vignette for details: https://github.com/mm-shah/MetaNeighbor/blob/master/vignettes/MetaNeighbor.pdf 
###
### Since Diaz-Mejia et al (2019) wanted to evaluate each cell cluster label assignment as 'testing' dataset,
### instead of their average as 'testing' and 'training', they followed advise from Maggie Crow (Metaneighbor developer)
### and Diaz-Mejia et al commented the following line in the MetaNeighbor source code `~/.../R/MetaNeighborUS.R`
### `cell_NV <- (cell_NV+t(cell_NV))/2`
### And then compiling MetaNeighbor from source in R with `install.packages()` as described below.
####################################

####################################
### Dependencies:
####################################
suppressPackageStartupMessages(library(optparse))             # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(SummarizedExperiment)) # (CRAN) to format inputs for SummarizedExperiment
suppressPackageStartupMessages(library(data.table))           # (CRAN) to speed up reading matrices with 'fread' instead of 'read.table'
suppressPackageStartupMessages(library(GSA))                  # (CRAN) to handle *gmt infile
### *** Package MetaNeighbor ***
### 1) Download the source from https://github.com/mm-shah/MetaNeighbor
### 2) Comment command `cell_NV <- (cell_NV+t(cell_NV))/2` in file `~/.../R/MetaNeighborUS.R`
### 3) Install package from source like: 
###    `install.packages("~/path_to/MetaNeighbor-master", repos = NULL, type="source")`
suppressPackageStartupMessages(library(MetaNeighbor))         
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
                ...etc
                Default = 'No default. It's mandatory to specify this parameter'"),
  
  make_option(c("-c", "--infile_signature"), default="NA",
              help="A path/name to a <tab> delimited *file* with the cell type signatures in either format:
                'gmt' format, like:
                GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 Gene3
                GeneSet2_ID  GeneSet2_Name  Gene4 Gene5
                ... etc

                OR

                'mat' gene expression profiles, like:
                GENE     GeneSet1_ID   GeneSet2_ID
                Gene1    0.91          0.0
                Gene2    0.92          0.0
                Gene3    0.93          0.0
                Gene4    0.0           0.94
                Gene5    0.0           0.95
                Note: values can also be binary (0's and 1's)
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-t", "--type_signature"), default="NA",
              help="Indicates if --infile_signature is either 'gmt' or 'mat' format
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-g", "--infile_cell_types"), default="NA",
              help="A path/name to a <tab> delimited *file* with the cell type annotations (positive gold standards), like:
              clust1  Cell_type_2
              clust2  Cell_type_3
              clust3  Cell_type_1
              Note: Cell_type ID's must match signature GeneSet ID's"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'")

)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileSignature <- opt$infile_signature
TypeSignature   <- opt$type_signature
InfileCellTypes <- opt$infile_cell_types
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles

Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################
writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infile_mat", "infile_signature", "type_signature", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

####################################
### Create outdirs
####################################
writeLines("\n*** Create outdirs ***\n")
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Tempdir)

### Checking if the system follows a mirror /scratch vs. /home structure
### If that's the case, this script will use Tempdir at /scratch, not at /home
### This is because systems like SciNet (Compute Canada) using 'Slurm Workload Manager'
### only allow to write in /scratch when using slave nodes, not in /home
### Comment the following '4' lines if you have a mirror structure and still want to use /home (called /Users in Mac)
### If you don't have a mirror structure just ignore this
ScratchTempdir<-gsub("home","scratch", Tempdir)
if (dir.exists(ScratchTempdir)[1] == 1) {
  Tempdir<-ScratchTempdir
}

dir.create(file.path(Outdir, "METANEIGHBOR"), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),        showWarnings = F, recursive = T)

####################################
### Load the matrix of genes (rows) vs. cell clusters (columns)
####################################
writeLines("\n*** Load the matrix of genes (rows) vs. cell clusters (columns) ***\n")

### Load the matrix of genes (rows) vs. cell clusters (columns)
InfileMat.df<-data.frame(fread(InfileMat, sep="\t", na.strings=c("NA")), row.names=1)

####################################
###  Load the matrix of genes (rows) vs. cell clusters (columns)
####################################
writeLines("\n*** Load the cell type signatures ***\n")

### Load the cell type signatures
if (length(grep('^gmt$', TypeSignature, perl = T, ignore.case = T))) {
  
  ### Creates a data.frame from gene sets
  gmt1<-GSA.read.gmt(InfileSignature)
  gmt2<-list()
  for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
  }
  AllGenesInSignature<-unique(unlist(x = gmt2, recursive = T, use.names = F))

  SigMat.df<-data.frame()
  for (gene_set_id in  names(gmt2)) {
    for (gene in AllGenesInSignature) {
      if (gene %in% gmt2[[gene_set_id]]) {
        SigMat.df[gene,gene_set_id] <- 1
      }else{
        SigMat.df[gene,gene_set_id] <- 0
      }
    }
  }

} else if (length(grep('^mat$', TypeSignature, perl = T, ignore.case = T))) {
  
  ### Creates a data.frame directly from signature matrix
  SigMat.df<-data.frame(fread(InfileSignature, sep="\t", na.strings=c("NA")), row.names=1)
  AllGenesInSignature<-rownames(SigMat.df)
  
} else {
  stop(paste("Parameter -c ", TypeSignature, " must be either  'gmt' or 'mat'. Use option -h for help.", sep = "", collapse = ""))
}
FullMat.df<-merge(SigMat.df, InfileMat.df, by = "row.names", all = T)
FullMat.df[is.na(FullMat.df)] <- 0
rownames(FullMat.df)<-FullMat.df[,"Row.names"]
FullMat.df$Row.names<-NULL

## Using study ID 'Signature' for all cell types from --infile_signature
## and study ID 'DataMatrix' for all cell clusters from --infile_mat
StudyIDs<-c(replicate(ncol(SigMat.df), "Signature"), replicate(ncol(InfileMat.df), "DataMatrix"))

####################################
### Load the cell type annotations
####################################
writeLines("\n*** Load the cell type annotations ***\n")

InfileCellTypes.df<-data.frame(fread(InfileCellTypes, sep="\t", header = F, select = c(1:2)), row.names = 1)
CellTypes<-c(colnames(SigMat.df), c(colnames(InfileMat.df)))

####################################
### Create a SummarizedExperiment object
####################################
writeLines("\n*** Create SummarizedExperiment object ***\n")

### Create SummarizedExperiment object
### 'chr', 'start', 'end' and 'strand' are columns needed to create the SummarizedExperiment object
FullMat.df$chr<-1
FullMat.df$start<-1
FullMat.df$end<-1
FullMat.df$strand<-"*"

FullMat.seo <- makeSummarizedExperimentFromDataFrame(FullMat.df)
FullMat.seo

####################################
### Run MetaNeighbor US
####################################
writeLines("\n*** Run MetaNeighbor US ***\n")

StartTimeMetaneighborUS<-Sys.time()
AUROC_scores_US = MetaNeighborUS(var_genes = AllGenesInSignature,
                             dat = FullMat.seo, 
                             study_id = StudyIDs,
                             cell_type = CellTypes)

EndTimeMetaneighborUS<-Sys.time()

### Print out METANEIGHBOR results
RowsDataMatrix<-grep(pattern = "^DataMatrix", rownames(AUROC_scores_US))
ColsDataMatrix<-grep(pattern = "^Signature",  colnames(AUROC_scores_US))
AUROC_scores_US_Datamatrix<-AUROC_scores_US[RowsDataMatrix,ColsDataMatrix]
colnames(AUROC_scores_US_Datamatrix)<-gsub(pattern = "^Signature\\|",  replacement = "", x = colnames(AUROC_scores_US_Datamatrix))
rownames(AUROC_scores_US_Datamatrix)<-gsub(pattern = "^DataMatrix\\|", replacement = "", x = rownames(AUROC_scores_US_Datamatrix))

OutfileEnrichmentScoresUS<-paste(Tempdir,"/",PrefixOutfiles,".MetaNeighborUS_AUROC.tsv", sep="")
Headers<-paste("MetaNeighborUS", paste(colnames(AUROC_scores_US_Datamatrix), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileEnrichmentScoresUS, row.names = F, col.names = F, sep="\t", quote = F)
write.table(AUROC_scores_US_Datamatrix, file = OutfileEnrichmentScoresUS, row.names = T, col.names = F, sep="\t", quote = F, append = T)

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".MetaNeighborUS_UsedOptions.txt", sep="")
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

TookTimeGsva    <-format(difftime(EndTimeMetaneighborUS, StartTimeMetaneighborUS, units = "secs"))
TookTimeOverall <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "secs"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".MetaNeighborUS_CPUusage.txt", sep="")
ReportTime<-c(
  paste("MetaNeighborUS",  TookTimeGsva,    collapse = "\t"),
  paste("overall",          TookTimeOverall, collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/METANEIGHBOR/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".MetaNeighborUS_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/METANEIGHBOR/",eachFile,sep=""),overwrite=T)
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

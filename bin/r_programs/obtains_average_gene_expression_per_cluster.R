####################################
### Script 'obtains_average_gene_expression_per_cluster.R' computes the average gene expression
### for each gene of --input, for each cluster of --input_clusters.
###
### --input can be either matrix market format files (barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz),
### or a sparse matrix with genes (rows) and cell IDs (columns)
###
### --input_clusters must be a <tab> delimited paired table with cell IDs and cluster IDs
###
### Will return a matrix with the average gene expression of each gene (rows) for each cell cluster (columns)
### 
### Follows this procedure, using the use.raw=T space https://satijalab.org/seurat/interaction_vignette.html
###
### NEW IMPLEMENTATIONS SINCE Seurat v2:
### 1) Rewritten with Seurat v3 commands (including ability to read output from Cell Ranger v3)
###    Main differences vs. Seurat v2 include:
###    a) new function names
###    b) since Cell Ranger v3 allows now to have multiple features (not only genes),
###       in general all references to 'genes' in v2 are now called 'features'
###
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_average_gene_expression_per_cluster.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
suppressPackageStartupMessages(library(Seurat))       # to run QC, differential gene expression and clustering analyses
suppressPackageStartupMessages(library(dplyr))        # needed by Seurat for data manupulation
suppressPackageStartupMessages(library(optparse))     # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(data.table))   # to read tables quicker than read.table - only needed is using '-t DGE'
####################################

####################################
### TO DO
### 1) In "Load data" we use Seurat library(Read10X). In this command, when all barcodes come from the same sample
### (i.e. barcode ID's finish with the same digit), like:
###    CTCTACGCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    CTGAAACCAAGAGGCT-1
###    ... etc
###    Read10X will remove the '-digit'
###
###    Hence we need to implement code to remove the '-digit' from --input_clusters barcode ID's as well
###    WHEN all barcodes come from the same sample.
###    For now, this script is removing the digit always. And user must provide the inputs like:
###    1-CTCTACGCAAGAGGCT
###    2-CTCGAAAAGCTAACAA
###    3-CTGCCTAGTGCAGGTA
###    Instead of:
###    CTCTACGCAAGAGGCT-1
###    CTCGAAAAGCTAACAA-2
###    CTGCCTAGTGCAGGTA-3
###    If their --input contains data from multiple samples
###   
####################################

####################################
### Turning warnings off for the sake of a cleaner aoutput
####################################
oldw <- getOption("warn")
options( warn = -1 )

####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-i", "--input"), default="NA",
              help="Either the path/name to a MTX *directory* with barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files;
                or path/name of a <tab> delimited digital gene expression (DGE) *file* with genes in rows vs. barcodes in columns
                Notes:
                The 'MTX' files can be for example the output from Cell Ranger 'count' v2 or v3: `/path_to/outs/filtered_feature_bc_matrix/`
                Cell Ranger v2 produces unzipped files and there is a genes.tsv instead of features.tsv.gz
                Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-t", "--input_type"), default="NA",
              help="Indicates if input is either a 'MTX' directory or a 'DGE' file
              Default = 'No default. It's mandatory to specify this parameter'"),
  #
  make_option(c("-c", "--input_clusters"), default="NA",
              help="A path/name to a <tab> delimited file, with one or more cluster labeling columns, in format like:
                Barcode            cluster_labels_1   cluster_labels_1 ... etc
                AAACCTGGTATTACCG   0                  0
                AAACCTGTCGGCATCG   1                  1
                AAACCTGTCTCCAACC   1                  2
                AAACCTGTCTTCGAGA   2                  3
                ... etc"),
  #
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  #
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID")
  )

opt <- parse_args(OptionParser(option_list=option_list))

Input          <- opt$input
InputType      <- opt$input_type
InputClusters  <- opt$input_clusters
Outdir         <- opt$outdir
PrefixOutfiles <- opt$prefix_outfiles

Tempdir        <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

####################################
### Define default parameters
####################################

DefaultParameters <- list(
  MinCells = 3,
  MinGenes = 200
)

####################################
### Start stopwatch
####################################

StartTimeOverall <-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("input", "input_type", "input_clusters", "outdir", "prefix_outfiles")
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
#
dir.create(file.path(Outdir, "AVERAGE_GENE_EXPRESSION"), recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

####################################
### Load scRNA-seq data
####################################
writeLines("\n*** Load scRNA-seq data ***\n")

if (regexpr("^MTX$", InputType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  input.matrix <- Read10X(data.dir = Input)
}else if (regexpr("^DGE$", InputType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  input.matrix <- data.frame(fread(Input),row.names=1)
}else{
  stop(paste("Unexpected type of input: ", InputType, "\n\nFor help type:\n\nRscript Runs_Seurat_Clustering.R -h\n\n", sep=""))
}
dim(input.matrix)

####################################
### Create a Seurat object
####################################
writeLines("\n*** Create a Seurat object ***\n")

seurat.object  <- CreateSeuratObject(counts = input.matrix, min.cells = DefaultParameters$MinCells, min.features = DefaultParameters$MinGenes, project = PrefixOutfiles)
seurat.object

####################################
### Load cluster assignments and get average expression
####################################
writeLines("\n*** Load cluster assignments and get average expression ***\n")

if (InputClusters == "NA") {
  stop("An infile provided by option -c was expected")
}else{
  CellClusters <- data.frame(read.table(InputClusters, header = T, row.names = 1, sep = "\t"))
  
  # This is because Seurat removes the last '-digit' from barcode ID's when all barcodes finish with the same digit (i.e. come from the same sample)
  # so that barcodes from --input_clusters and --input can match each other
  rownames(CellClusters)<-gsub(x =rownames(CellClusters), pattern = "-[0-9]+$", perl = T, replacement = "")
  seurat.object <- AddMetaData(object = seurat.object, metadata = CellClusters)
}

for (cluste_type in colnames(CellClusters)) {
  
  print(paste("Getting expression averages for: ", cluste_type, sep=""))
  
  # switch the identity class of all cells to reflect "clusters"
  Idents(object = seurat.object) <- cluste_type

  ####################################
  ### Get average expression for each cluster for each gene
  ####################################
  cluster.averages<-AverageExpression(object = seurat.object, use.counts = T)

  OutfileClusterAverages<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_", cluste_type, ".tsv", sep="")
  
  if ((grepl(pattern = "^[0-9]+$", x = names(cluster.averages)[[1]])) == TRUE) {
    Headers<-paste("AVERAGE_GENE_EXPRESSION",paste("C", names(cluster.averages[["RNA"]]), sep="",collapse="\t"),sep="\t",collapse = "\t")
  }else{
    Headers<-paste("AVERAGE_GENE_EXPRESSION",paste(names(cluster.averages[["RNA"]]),sep="",collapse="\t"),sep="\t",collapse = "\t")
  }
  write.table(Headers,file = OutfileClusterAverages, row.names = F, col.names = F, sep="\t", quote = F)
  write.table(data.frame(cluster.averages),file = OutfileClusterAverages, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_UsedOptions.txt", sep="")
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

TookTimeOverall        <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "min"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".SEURAT_AverageGeneExpressionPerCluster_CPUusage.txt", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"/AVERAGE_GENE_EXPRESSION/",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".SEURAT_AverageGeneExpressionPerCluster", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir,"/AVERAGE_GENE_EXPRESSION/",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"/",eachFile,sep=""))
})

print(paste("Outfiles at: ", Outdir, "/AVERAGE_GENE_EXPRESSION/", sep = "", collapse = ""))

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



####################################
### Script 'propagates_permuted_gmt_files_to_profile.R' takes as inputs:
### a) the *PERMUTE_GMT_UsedOptions.txt log file from script 'obtains_permuted_samples_from_gmt.R'
### b) a --infile_profile with gene expression profiles
### and propagates the subsampled gene sets from 'obtains_permuted_samples_from_gmt.R'
### into new gene expression profiles, where genes that were removed from subsampled gene sets
### are replaced in --infile_profile by a --value_replacing (e.g. min, median, mean or 0.0)
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/propagates_permuted_gmt_files_to_profile.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
### 'optparse'   (CRAN) to handle one-line-commands
### 'GSA'        (CRAN) to handle *gmt infile
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GSA))
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
  make_option(c("-c", "--infile_profile"), default="NA",
              help="A path/name to a <tab> delimited gene expression profile matrix, like:
              GENE   CellType1  CellType2  CellType3 ... etc
              Gene1  0.1        0.04       0.6
              Gene2  0.02       0.1        0.01
              Gene3  0.04       0.3        0.06
              ...etc"),
  
  make_option(c("-l", "--infile_log_permutted_gmt"), default="NA",
              help="A path/name to a <tab> delimited *PERMUTE_GMT_UsedOptions.txt log file from 'obtains_permuted_samples_from_gmt.R'
              This will be used to get: --infile_gmt, --iterations, --sample_percentage, and --prefix_outfiles options"),
  
  make_option(c("-v", "--value_replacing"), default="median",
              help="Script 'obtains_permuted_samples_from_gmt.R' removes genes from each cell-type (gene-set)
              Whereas here, rows (genes) and columns (cell-types) from --infile_profile will remain the same,
              but values of genes removed by 'obtains_permuted_samples_from_gmt.R' will be replaced by --value_replacing
              You can use a custom --value_replacing (e.g. '0.0' or 'NA'), or type 'median', 'mean' or 'min'
              to compute and use such metric for replacements on each column of --infile_profile"),
  
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileProfile          <- opt$infile_profile
InfileLogPermuttedGmt  <- opt$infile_log_permutted_gmt
ValueReplacing         <- opt$value_replacing
Outdir                 <- opt$outdir

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_profile", "infile_log_permutted_gmt", "value_replacing", "outdir")
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
Outdir<-gsub("/$", "", Outdir)
#
dir.create(file.path(Outdir, "PROFILE_PERMUTATIONS"), showWarnings = F, recursive = T)

####################################
### Get options from *PERMUTE_GMT_UsedOptions.txt log file
####################################
writeLines("\n*** Get options from *PERMUTE_GMT_UsedOptions.txt log file ***\n")

lines <-readLines(file(InfileLogPermuttedGmt, open="r"))
for (i in 1:length(lines)) {
  fields = strsplit(lines[i], '\t') [[1]]
  
  if (grepl(pattern = "^-c", ignore.case = T, x = lines[i]) == T) {
    InfileFullGmt     <- fields[3]
  }else if (grepl(pattern = "^-r", ignore.case = T, x = lines[i]) == T) {
    Iterations        <- fields[3]
  }else if (grepl(pattern = "^-s", ignore.case = T, x = lines[i]) == T) {
    SamplePercentages <- fields[3]
  }else if (grepl(pattern = "^-p", ignore.case = T, x = lines[i]) == T) {
    PrefixOutfiles    <- fields[3]
  }else if (grepl(pattern = "^-o", ignore.case = T, x = lines[i]) == T) {
    OutdirReducedGmt  <- fields[3]
  }
}
close(file(InfileLogPermuttedGmt))

####################################
### Check that --sample_percentages and --iterations options are defined correctly
####################################

if (grepl(pattern = "^[0-9]+,[0-9]+,[0-9]+$", ignore.case = T, x = SamplePercentages) == T) {
  ListPercentLimits<-unlist(strsplit(SamplePercentages, ","))
  MinPercent<-as.numeric(ListPercentLimits[1])
  MaxPercent<-as.numeric(ListPercentLimits[2])
  PercentBy <-as.numeric(ListPercentLimits[3])
}else{
  stop(paste("Unexpected format of --sample_percentages '-s ",  SamplePercentages, "'. Expected format is like '-s 10,100,10'", sep = "", collapse = ""))
}

if (grepl(pattern = "^[0-9]+-[0-9]+$", ignore.case = T, x = Iterations) == T) {
  IterationLimits<-unlist(strsplit(Iterations, "-"))
  MinIterations<-as.numeric(IterationLimits[1])
  MaxIterations<-as.numeric(IterationLimits[2])
}else{
  stop(paste("Unexpected format of --iterations '-r ",  Iterations, "'. Expected format is like '-r 1-100'", sep = "", collapse = ""))
}

####################################
### Get values to replace
####################################
writeLines("\n*** Get values to replace ***\n")

OriginalMat<-read.table(InfileProfile,header = T, row.names = 1, sep = "\t")

### Get values to replace
ValuesReplacing<-list()
if (grepl(pattern = "^median$|^mean$|^min$", ignore.case = T, x = ValueReplacing) == T) {
  ValuesReplacing<-apply(X = OriginalMat, MARGIN = 2, FUN = ValueReplacing)
}else{
  for (colName in colnames(OriginalMat)) {
    ValuesReplacing[[colName]] <- ValueReplacing
  }
}

####################################
###  Load full gmt file
####################################
writeLines("\n*** Load full gmt file ***\n")

### Creates FullGmt object with the FULL gene set memberships
gmt1<-GSA.read.gmt(InfileFullGmt)
FullGmt<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  FullGmt[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

####################################
###  Load reduced gmt files and outputs new profiles
####################################
writeLines("\n*** Load reduced gmt files and outputs new profiles ***\n")

for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
  print(paste("Generate reduced files for percent:", percent, sep = " ", collapse = " "))
  for (perm in c(MinIterations:MaxIterations)) {

    ### Load reduced gmt file
    InfileReducedGmt<-paste(OutdirReducedGmt, "/GMT_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="")
    gmt2<-GSA.read.gmt(InfileReducedGmt)
    ReducedGmt<-list()
    for (i in 1:length(gmt2[[1]])){
      tmp<-unlist(gmt2[[1]][i])
      ReducedGmt[[gmt2[[3]][i]]]<-tmp[which(tmp!="")]
    }
    
    ### Get genes to replace values per column
    NewMat<-OriginalMat
    for (colName in colnames(OriginalMat)) {
      GenesInFullGmt    <- FullGmt[[colName]]
      GenesInReducedGmt <- ReducedGmt[[colName]]
      GenesToReplaceValues <- setdiff(GenesInFullGmt,GenesInReducedGmt)
      NewMat[c(GenesToReplaceValues),colName]<-ValuesReplacing[[colName]]
    }
    
    ### Prepare data for permutted outfile
    OutfileNewMat<-paste(Outdir,"/PROFILE_PERMUTATIONS/",PrefixOutfiles, ".percent", percent, ".perm", perm, ".profile.tsv", sep="")
    file.create(file = OutfileNewMat)
    
    Headers<-paste("GENE",paste(colnames(NewMat),sep="",collapse="\t"), sep="\t", collapse = "\t")
    write.table(Headers,file = OutfileNewMat, row.names = F, col.names = F, sep="\t", quote = F)
    write.table(NewMat, file = OutfileNewMat, row.names = T, col.names = F, sep="\t", quote = F, append = T)
  }
}

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_PROFILE_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Report time used
####################################
EndTimeOverall<-Sys.time()

TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "secs"))

OutfileCPUusage<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_PROFILE_CPUusage.txt", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

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

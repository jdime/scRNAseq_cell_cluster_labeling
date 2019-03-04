####################################
### Script 'obtains_permuted_samples_from_gmt.R' subsamples each gene set from --infile_gmt
### retaining percentages of genes indicated by --sample_percentages
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_permuted_samples_from_gmt.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
suppressPackageStartupMessages(library(optparse)) # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(GSA))      # (CRAN) to handle *gmt infile
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
  make_option(c("-c", "--infile_gmt"), default="NA",
              help="A path/name to a <tab> delimited *file* of gene sets in *gmt format, like:
              GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 ... etc
              GeneSet2_ID  GeneSet2_Name  Gene2 Gene3 ... etc"),
  
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  
  make_option(c("-r", "--iterations"), default="1-100",
              help="Indicates the number of times to subsample --infile_gmt, e.g. '1-100' to run from iteration 1 to 100,
              or '50-100' to run from iteration 50 to 100
              Note: if using '-t N', then this script will use previously permuted *gmt infiles"),
  
  make_option(c("-s", "--sample_percentages"), default="10,50,10",
              help="Indicates three values: 'minimum_percentage', 'maximum_percentage' and 'increment' to randomly select a percent of genes from each gene set.
              For example, to select from 10% up to 50%, increasing by 5% (10%, 15%, 20%, ... 50%) use '10,50,5'"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              Note this script will automatically add 'percentXX.permYY' to the outfile name indicating replaced percentage and iteration number")
  
  )

opt <- parse_args(OptionParser(option_list=option_list))

InfileGmt         <- opt$infile_gmt
Outdir            <- opt$outdir
PrefixOutfiles    <- opt$prefix_outfiles
Iterations        <- opt$iterations
SamplePercentages <- opt$sample_percentage

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_gmt", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

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
### Create outdirs and define outfiles
####################################
writeLines("\n*** Create outdirs ***\n")
CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Outdir<-gsub("/$", "", Outdir)
#
dir.create(file.path(Outdir, "GMT_PERMUTATIONS"), showWarnings = F, recursive = T)

####################################
### Get GMT file permutations
####################################
writeLines("\n*** Get GMT file permutations ***\n")

### Creates OriginalGmt object with the gene set memberships
gmt1<-GSA.read.gmt(InfileGmt)
OriginalGmt<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  OriginalGmt[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
  for (perm in c(MinIterations:MaxIterations)) {
    
    ### Open permuted gmt outfile
    OutfileName<-paste(Outdir,"/GMT_PERMUTATIONS/",PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="")
    file.create(file = OutfileName)
    
    classNumber<-0
    for (class in OriginalGmt) {
      classNumber<-classNumber+1
      NtoGet<-length(class) * percent / 100
      if (NtoGet < 1) {
        NtoGet<-1
      }
      newClass<-paste(c(gmt1$geneset.names[classNumber], gmt1$geneset.descriptions[classNumber],sample(x=OriginalGmt[[classNumber]], size = NtoGet)),sep="\t",collapse = "\t")
      write.table(newClass,file = OutfileName, row.names = F, col.names = F, sep="\t", quote = F, append = T)
    }
  }
}

####################################
### Report used options
####################################
OutfileOptionsUsed<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_GMT_UsedOptions.txt", sep="")
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

OutfileCPUusage<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_GMT_CPUusage.txt", sep="")
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

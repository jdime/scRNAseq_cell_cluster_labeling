####################################
### Script 'subsamples_gene_classes_and_runs_enrichment_scripts.R' is a *wrapper* to run implementations of
### methods to predict cell types using cell clusters from scRNA-seq data and cell type gene expression signatures
###
### It computes:
### a) generating Receiver Operating Characteristic (ROC) and Precision-Recall (PR) curves for each method predictions
### b) getting ROC and PR Area Under the Curves (ROC AUC, and PR AUC)
### c) running robustness analyses by subsampling marker genes from cell type signatures and repeating ROC AUC and PR AUC analyses
###
### Uses three types of infiles:
### (1) a matrix with average gene expression for each gene (rows) cell clusters from scRNA-seq (columns)
### (2) cell type gene expression signatures
### (3) a reference paired table with gold standard cell type predictions for each cell cluster
###
### The general flow of this script is the following:
### (i)   Formats inputs for each requested cell type prediction method: CIBERSORT, GSEA, GSVA, METANEIGHBOR and ORA
### (ii)  Subsamples cell type gene expression signatures using infiles (2)
### (iii) Runs cell type predictions methods using infile (1) and subsampled signatures from step (ii)
### (iv)  Gather results from step (iii)
### (v)   Runs ROC, PR and AUC analyses using results from step (iv) and infile (3)
### (vi)  Runs robustness analyses (optional), including violin plots of AUC distributions with results from step (v)
###
### Notes:
###        Note 1
###        In the cases of CIBERSORT and METANEIGHBOR, two forms of cell type signature infiles can be provided:
###        - in the form of gene sets (2a) or in the form of gene expression profiles (2b)
###        - Subsampling genes from is always conducted using 2a type files, which can be propagated to 2b type files using script:
###          'propagates_permuted_gmt_files_to_profile.R'
###
###        Note 2
###        To add new methods to evaluate you can:
###        a) add the path to the executable in 'Dependencies:'
###        b) provide inputs and parameters adding commands indicated by 'Runs NEW_METHOD and gather its outputs'
###        c) add program name to this script 'option_list' parameters --software_to_run and --software_to_auc
###        d) add program name to script 'obtains_performance_plots_from_cluster_labelings.pl' so that it can load its predictions outfile
###        e) assign a colour from list(ColourDefinitions) to the new method at list(ListColoursViolin)
###
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/subsamples_gene_classes_and_runs_enrichment_scripts.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
### R and Rscript (tested with version 3.5.1)
###
### R libraries:
suppressPackageStartupMessages(library(optparse)) # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(vioplot))  # (CRAN) to generate violin plots of permutation results
###
### External scripts (check each script for their own dependencies)
### and change their path as needed here:
PathNameToPermuteGmtRscript       <- "~/r_programs/obtains_permuted_samples_from_gmt.R"
PathNameToPermuteProfileRscript   <- "~/r_programs/propagates_permuted_gmt_files_to_profile.R"
PathNameToCibersortscript         <- "~/bin/obtains_CIBERSORT_for_MatrixColumns.pl"
PathNameToGSEAscript              <- "~/bin/obtains_GSEA_for_MatrixColumns.pl"
PathNameToGSVAscript              <- "~/r_programs/obtains_GSVA_for_MatrixColumns.R"
PathNameToMETANEIGHBORscript      <- "~/r_programs/obtains_METANEIGHBOR_for_MatrixColumns.R"
PathNameToORAscript               <- "~/bin/obtains_ORA_for_MatrixColumns.pl"
PathNameToPerformancePlotsAndAucs <- "~/bin/obtains_performance_plots_from_cluster_labelings.pl"
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
              help="Path/name to a <tab> delimited *file* with average gene expression per cell cluster from scRNA-seq,
                genes in rows and cell clusters in columns, like:
                genes  clust1  clust2  clust3 ... etc
                Gene1  0       0.0045  0.0008
                Gene2  0.0077  0.0175  0.0082
                Gene3  0.0800  0.1532  0.0745
                ...etc
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-c", "--infile_signature_gmt"), default="NA",
              help="Path/name to a <tab> delimited cell type signature *file* in the form of gene lists (*gmt format), like:
                GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 ... etc
                GeneSet2_ID  GeneSet2_Name  Gene2 Gene3 ... etc
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-g", "--infile_gold"), default="NA",
              help="A path/name to a <tab> delimited *file* of gold standard cluster labels in format, like:
                clust1  GeneSet1_ID
                clust2  GeneSet3_ID
                clust3  GeneSet5_ID
                ... etc
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-t", "--permute_gmt"), default="n",
              help="Indicates if permutation of --infile_signature_gmt should be conducted [use y/Y] or permuted files exist already [type n/N]
                Default = 'n'"),
  
  make_option(c("-u", "--propagate_permuted_gmt_to_profiles"), default="n",
              help="Only needed by 'CIBERSORT_PROFILE' and 'METANEIGHBOR_PROFILE'
                Indicates if permutation of --infile_signature_gmt files should be propagated to profiles [use y/Y] or permuted files exist already [type n/N]
                Default = 'n'"),

  make_option(c("-d", "--infile_signature_profile"), default="NA",
              help="Only needed by 'CIBERSORT_PROFILE' and 'METANEIGHBOR_PROFILE'
                Path/name to a <tab> delimited cell type signature *file* in the form of gene expression profile matrix, like:
                GENE   CellType1  CellType2  CellType3 ... etc
                Gene1  0.1        0.04       0.6
                Gene2  0.02       0.1        0.01
                Gene3  0.04       0.3        0.06
                ...etc
                Or type 'NA' if neither 'CIBERSORT_PROFILE' nor 'METANEIGHBOR_PROFILE' methods are being used"),
              
  make_option(c("-s", "--sample_percentages"), default="10,100,10",
              help="Indicates three values: minimum_percentage, maximum_percentage, increment to be used to sample --infile_signature_gmt
                For example, to subsample percentages from 10 up to 50, increasing by 5 (i. e. 10, 15, 20, ... 50) use '10,50,5'
                Note: if using '-n N' then this script will take previously permuted *gmt infiles, stored at /path_to/--outdir/GMT_PERMUTATIONS
                Default = '10,100,10'"),
  
  make_option(c("-r", "--iterations"), default="1-100",
              help="Indicates the number of times to subsample --infile_signature_gmt
                For example, to run iterations 1 to 100 use '1-100' or to run iteration 50 to 100 use '50-100'
                Note: if using '-t N', then this script will use previously permuted *gmt infiles
                Default = '1-100'"),
  
  make_option(c("-m", "--software_to_run"), default="NONE",
              help="Indicates <comma> delimited name(s) of enrichment sofware to *run*:
                CIBERSORT_BINARY,CIBERSORT_PROFILE,GSEA,GSVA,METANEIGHBOR_BINARY,METANEIGHBOR_PROFILE,ORA
                Default = 'NONE'."),
  
  make_option(c("-n", "--software_to_auc"), default="NONE",
              help="Indicates <comma> delimited name(s) of enrichment sofware results to *gather* together:
                CIBERSORT_BINARY,CIBERSORT_PROFILE,GSEA,GSVA,METANEIGHBOR_BINARY,METANEIGHBOR_PROFILE,ORA
                This may be useful to gather results from separate runs
                Default = 'NONE'."),
  
  make_option(c("-v", "--generate_violin_plots"), default="n",
              help="Indicates if violin plots should be generated from --software_to_run and --software_to_auc steps. Type [y/Y] or [n/N]
                Default = 'n'"),
  
  make_option(c("-a", "--roc_auc_violin_plot_y_axes_limits"), default="0.5,1",
              help="Indicates the min,max values for the ROC AUC violin plot y-axes of enrichment software
                Default = '0.5,1'"),
  
  make_option(c("-b", "--pr_auc_violin_plot_y_axes_limits"), default="0,1",
              help="Indicates the min,max values for the PR AUC violin plot y-axes of enrichment software
                Default = '0,1'"),
  
  make_option(c("-k", "--roc_horizontal_guide_line"), default="NA",
              help="Indicates a value in the y-axis of the ROC AUC violin plots to place a horizontal guide line. Or type 'NA' to omit
                Default = 'NA'"),
  
  make_option(c("-l", "--pr_horizontal_guide_line"), default="NA",
              help="Indicates a value in the y-axis of the PR AUC violin plots to place a horizontal guide line. Or type 'NA' to omit
                Default = 'NA'"),
  
  make_option(c("-w", "--print_plot_ticks_labels_software"), default="all",
              help="Indicates if plots should contain:
                a) tick marks, axes values, axes labels and enrichment sofware name (type 'all')
                b) only tick marks (type 'ticks')
                c) only tick marks and axes values (type 'tick_values')
                d) only tick marks, axes values and axes labels (type 'tick_val_lab')
                e) only enrichment sofware name (type 'software')
                f) no tick marks, axes values, axes labels and enrichment sofware name (type 'none')
                Default = 'all'"),

  make_option(c("-e", "--plot_height"), default="1",
              help="Indicates a factor to decrease/increase the height of the plot
                For example, to use a height 1.5 times bigger than normal use '-e 1.5', or to reduce it by half use '-e 0.5'
                Default = '1'"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Note this script will automatically add 'percentXX.permYY' to the outfile name indicating -s value and iteration number
                Default = 'No default. It's mandatory to specify this parameter'"),
  
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-y", "--ora_use_values_or_ranks"), default="NA",
              help="This applies only to the ORA test:
                Indicates if a cutoff applied to --infile_mat to get 'white/drawn' balls should be based on scores 'ranks' or 'values'
                For example, using '-y values -z 1.5' will use all genes with values >= 1.5 in each column of --infile_mat as 'drawn' balls
                Whereas, using '-y ranks -z 100' will sort each column of --infile_mat and use the top 100 ranked genes as 'drawn' balls
                Or type 'NA' if ORA method is not being used
                Default = 'NA'"),
  
  make_option(c("-z", "--ora_mat_cutoff"), default="NA",
              help="This applies only to the ORA test:
                Indicates the cutoff to be applied to --infile_mat to get 'white/drawn' balls
                Or type 'NA' if ORA method is not being used
                Default = 'NA'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat            <- opt$infile_mat
InfileGmt            <- opt$infile_signature_gmt
InfileGold           <- opt$infile_gold
PermuteGmt           <- opt$permute_gmt
PermuteProfile       <- opt$propagate_permuted_gmt_to_profiles
InfileProfile        <- opt$infile_signature_profile
SamplePercentages    <- opt$sample_percentages
Iterations           <- opt$iterations
SoftwareToRun        <- opt$software_to_run
SoftwareToAuc        <- opt$software_to_auc
PlotViolin           <- opt$generate_violin_plots
AverageAucsTable     <- opt$generate_average_auc_tables
RocAucAxes           <- opt$roc_auc_violin_plot_y_axes_limits
PrAucAxes            <- opt$pr_auc_violin_plot_y_axes_limits
PrintLabelsAndTicks  <- opt$print_plot_ticks_labels_software
PlotHeight           <- as.numeric(opt$plot_height)
RocAucAbbline        <- opt$roc_horizontal_guide_line
PrAucAbbline         <- opt$pr_horizontal_guide_line
PrefixOutfiles       <- opt$prefix_outfiles
Outdir               <- opt$outdir
OraValOrRank         <- opt$ora_use_values_or_ranks
OraCutoff            <- opt$ora_mat_cutoff

StartTimeOverall<-Sys.time()

####################################
### Define default parameters
####################################

### Colour definitions
ColourDefinitions<-list("orange"        = "E69F00",
                        "bluishgreen"   = "009E73",
                        "reddishpurple" = "CC79A7",
                        "dodgerblue"    = "1E90FF",
                        "vermillion"    = "D55E00",
                        "snow4"         = "8B8989",
                        "yellow"        = "FFD700",
                        "seagreen3"     = "43CD80",
                        "pink"          = "FFC0CB",
                        "skyblue"       = "56B4E9",
                        "orchid4"       = "8B4789",
                        "blue"          = "0072B2",
                        "black"         = "000000"
)

ListColoursViolin<-list("CIBERSORT_BINARY"     = ColourDefinitions[["orange"]][[1]],
                        "CIBERSORT_PROFILE"    = ColourDefinitions[["bluishgreen"]][[1]],
                        "GSEA"                 = ColourDefinitions[["reddishpurple"]][[1]],
                        "GSVA"                 = ColourDefinitions[["dodgerblue"]][[1]],
                        "METANEIGHBOR_BINARY"  = ColourDefinitions[["yellow"]][[1]],
                        "METANEIGHBOR_PROFILE" = ColourDefinitions[["pink"]][[1]],
                        "ORA"                  = ColourDefinitions[["vermillion"]][[1]]
)

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_mat", "sample_percentages", "iterations", "outdir", "prefix_outfiles")
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
### Create outdirs and define outfiles and global variables
####################################
writeLines("\n*** Create outdirs ***\n")

CommandsToGetUserHomeDirectory<-("eval echo \"~$USER\"")
UserHomeDirectory<-system(command = CommandsToGetUserHomeDirectory, input = NULL, wait = T, intern = T)
#
Outdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Outdir<-gsub("/$", "", Outdir)

### This file will contain the list of outfiles to be used to get AUCs
ListExpectedOutfiles<-list()
OutfileListToGetAucs<-paste(Outdir, "/", PrefixOutfiles, ".PERMUTE.listToGetAuc.tsv", sep="", collapse = "")
file.create(OutfileListToGetAucs)

####################################
### Define external script paths
####################################
writeLines("\n*** Search for external dependencies ***\n")

ListExternalScripts<-list("PathNameToPermuteGmtRscript"       = PathNameToPermuteGmtRscript,
                          "PathNameToPermuteProfileRscript"   = PathNameToPermuteProfileRscript,
                          "PathNameToCibersortscript"         = PathNameToCibersortscript,
                          "PathNameToGSEAscript"              = PathNameToGSEAscript,
                          "PathNameToGSVAscript"              = PathNameToGSVAscript,
                          "PathNameToMETANEIGHBORscript"      = PathNameToMETANEIGHBORscript,
                          "PathNameToORAscript"               = PathNameToORAscript,
                          "PathNameToPerformancePlotsAndAucs" = PathNameToPerformancePlotsAndAucs
)

### Looking that external scripts exist
for (s in names(ListExternalScripts)) {
  ListExternalScripts[s]<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), ListExternalScripts[s])
  if (file.exists(ListExternalScripts[[s]])) {
    print(paste("Ok found: ", ListExternalScripts[[s]], sep = "", collapse = ""))
  } else {
    stop(paste("Couldn't find ", ListExternalScripts[[s]], " . Check section 'Define external script paths'", sep="", collapse = ""))
  }
}

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

####################################
### Get GMT file permutations
####################################
writeLines("\n*** Get GMT file permutations ***\n")

if (grepl(pattern = "y", ignore.case = T, x = PermuteGmt) == T) {
  CommandsToGetIterations<-paste("Rscript ", ListExternalScripts["PathNameToPermuteGmtRscript"],
                                 " -c ", InfileGmt,
                                 " -o ", Outdir,
                                 " -r ", Iterations,
                                 " -s ", SamplePercentages,
                                 " -p ", PrefixOutfiles,
                                 sep = "", collapse = "")
  
  system(CommandsToGetIterations, wait = T)
}else{
  print(c("Not permutting --infile_signature_gmt. Assuming permuted *gmt files exist at: ", paste(Outdir, "/GMT_PERMUTATIONS/", sep = "", collapse = "")))
}

####################################
### Propagate GMT permutations to Profiles
####################################

if (grepl(pattern = "y", ignore.case = T, x = PermuteProfile) == T) {
  writeLines("\n*** Propagate GMT permutations to Profiles ***\n")

  InfileGmtPermutationLog<-paste(Outdir, "/", PrefixOutfiles, ".PERMUTE_GMT_UsedOptions.txt", sep="", collapse = "")
  if(file.exists(InfileGmtPermutationLog)) {
    CommandsToRunPropagation<-paste("Rscript ", ListExternalScripts["PathNameToPermuteProfileRscript"],
                                    " -c ", InfileProfile,
                                    " -l ", InfileGmtPermutationLog,
                                    " -v ", "min",
                                    " -o ", Outdir,
                                    sep = "", collapse = "")
        
  ### Runs GMT propagations to Profiles
  system(CommandsToRunPropagation, wait = T)

  }else{
    stop(paste("Couldn't find GMT propagation log file: '", InfileGmtPermutationLog,  "'",sep="", collapse = ""))
  }
}else if (grepl(pattern = "_PROFILE", ignore.case = T, x = SoftwareToRun) == T) {
  print(c("Not permutting --infile_signature_profile. Assuming permuted gene expression profile files already exist at: ", paste(Outdir, "/PROFILE_PERMUTATIONS/", sep = "", collapse = "")))
}

####################################
### Runs CIBERSORT_BINARY and gather its outputs
####################################

ProgramName       <- "CIBERSORT_BINARY"
UsesGmtOrProfile  <- "GMT"
GmtOrProfileSuffix  <- ".gmt"
OutdirFromScript  <- "CIBERSORT"
SuffixPredictions <- ".CIBERSORT_enrichment_scores.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunCIBERSORT_BINARY<-paste(ListExternalScripts["PathNameToCibersortscript"],
                                           " -infile_matrix ",   InfileMat,
                                           " -infile_classes ",  InfileClassesLocal,
                                           " -path_outfiles ",   Outdir, "/", ProgramName, "_PERMUTATIONS/",
                                           " -prefix_outfiles ", PrefixOutfileLocal,
                                           " -nperm ",           1000,
                                           " -classes_type gmt ",
                                           sep = "", collapse = "")
      
      ### Runs CIBERSORT_BINARY
      system(CommandsToRunCIBERSORT_BINARY, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs CIBERSORT_PROFILE and gather its outputs
####################################

ProgramName       <- "CIBERSORT_PROFILE"
UsesGmtOrProfile  <- "PROFILE"
GmtOrProfileSuffix  <- ".profile.tsv"
OutdirFromScript  <- "CIBERSORT"
SuffixPredictions <- ".CIBERSORT_enrichment_scores.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunCIBERSORT_PROFILE<-paste(ListExternalScripts["PathNameToCibersortscript"],
                                            " -infile_matrix ",   InfileMat,
                                            " -infile_classes ",  InfileClassesLocal,
                                            " -path_outfiles ",   Outdir, "/", ProgramName, "_PERMUTATIONS/",
                                            " -prefix_outfiles ", PrefixOutfileLocal,
                                            " -nperm ",           1000,
                                            " -classes_type profile ",
                                            sep = "", collapse = "")
      
      ### Runs CIBERSORT_PROFILE
      system(CommandsToRunCIBERSORT_PROFILE, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")

  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs GSEA and gather its outputs
####################################

ProgramName         <- "GSEA"
UsesGmtOrProfile    <- "GMT"
GmtOrProfileSuffix  <- ".gmt"
OutdirFromScript    <- "GSEA"
SuffixPredictions   <- ".GSEA.Pval.Unfiltered.pos.mLog.mat.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunGSEA<-paste(ListExternalScripts["PathNameToGSEAscript"],
                               " -infile_matrix ",   InfileMat,
                               " -infile_classes ",  InfileClassesLocal,
                               " -path_outfiles ",   Outdir, "/", ProgramName, "_PERMUTATIONS/",
                               " -prefix_outfiles ", PrefixOutfileLocal,
                               " -cutoff_print_p ",  0.01,
                               " -cutoff_print_q ",  0.25,
                               " -nperm ",           1000,
                               " -use_universe i",
                               sep = "", collapse = "")
      
      ### Runs GSEA
      system(CommandsToRunGSEA, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs GSVA and gather its outputs
####################################
ProgramName         <- "GSVA"
UsesGmtOrProfile    <- "GMT"
GmtOrProfileSuffix  <- ".gmt"
OutdirFromScript    <- "GSVA"
SuffixPredictions   <- ".GSVA_enrichment_scores.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunGSVA<-paste("Rscript ", ListExternalScripts["PathNameToGSVAscript"],
                               " -i ", InfileMat,
                               " -c ", InfileClassesLocal,
                               " -o ", Outdir, "/", ProgramName, "_PERMUTATIONS/",
                               " -p ", PrefixOutfileLocal,
                               " -e 0.05 -f 0.1",
                               sep = "", collapse = "")
      ### Runs GSVA
      system(CommandsToRunGSVA, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs METANEIGHBOR_BINARY and gather its outputs
####################################

ProgramName       <- "METANEIGHBOR_BINARY"
UsesGmtOrProfile  <- "GMT"
GmtOrProfileSuffix  <- ".gmt"
OutdirFromScript  <- "METANEIGHBOR"
SuffixPredictions <- ".MetaNeighborUS_AUROC.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunMETANEIGHBOR_BINARY<-paste("Rscript ", ListExternalScripts["PathNameToMETANEIGHBORscript"],
                                           " -i ",  InfileMat,
                                           " -c ",  InfileClassesLocal,
                                           " -g ",  InfileGold,
                                           " -o ",  Outdir, "/", ProgramName, "_PERMUTATIONS/",
                                           " -p ",  PrefixOutfileLocal,
                                           " -t gmt ",
                                           sep = "", collapse = "")
      
      ### Runs METANEIGHBOR_BINARY
      system(CommandsToRunMETANEIGHBOR_BINARY, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs METANEIGHBOR_PROFILE and gather its outputs
####################################

ProgramName         <- "METANEIGHBOR_PROFILE"
UsesGmtOrProfile    <- "PROFILE"
GmtOrProfileSuffix  <- ".profile.tsv"
OutdirFromScript    <- "METANEIGHBOR"
SuffixPredictions   <- ".MetaNeighborUS_AUROC.tsv"

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")

      CommandsToRunMETANEIGHBOR_PROFILE<-paste("Rscript ", ListExternalScripts["PathNameToMETANEIGHBORscript"],
                                               " -i ",  InfileMat,
                                               " -c ",  InfileClassesLocal,
                                               " -g ",  InfileGold,
                                               " -o ",  Outdir, "/", ProgramName, "_PERMUTATIONS/",
                                               " -p ",  PrefixOutfileLocal,
                                               " -t mat ",
                                               sep = "", collapse = "")
      
      ### Runs METANEIGHBOR_PROFILE
      system(CommandsToRunMETANEIGHBOR_PROFILE, wait = T)
    }
  }
}


if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, SuffixPredictions, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs ORA and gather its outputs
####################################
ProgramName       <- "ORA"
UsesGmtOrProfile  <- "GMT"
GmtOrProfileSuffix  <- ".gmt"
OutdirFromScript  <- "ORA"
SuffixPredictions <- paste(".ORA.Pval.Uncorrected.cutoff_pos.mLog.mat.tsv", sep = "", collapse = "")

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToRun) == T) {
  cat("\n*** Running ", ProgramName,  " ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Running ", ProgramName, " : with ", percent, "% of data - iteration ", perm, "\n")
      #
      InfileClassesLocal<-paste(Outdir,"/", UsesGmtOrProfile, "_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, GmtOrProfileSuffix, sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")

      CommandsToRunORA<-paste(ListExternalScripts["PathNameToORAscript"],
                                 " -infile_matrix ",       InfileMat,
                                 " -infile_classes ",      InfileClassesLocal,
                                 " -path_outfiles ",       Outdir, "/", ProgramName, "_PERMUTATIONS/",
                                 " -prefix_outfiles ",     PrefixOutfileLocal,
                                 " -cutoff_pos ",          OraCutoff,
                                 " -use_values_or_ranks ", OraValOrRank,
                                 " -cutoff_ora ",          "0.05",
                                 " -cutoff_neg  ",         "NA",
                                 " -p_correction_test ",   "BH",
                                 " -use_ora ",             "pc",
                                 " -use_universe ",        "i",
                                 sep = "", collapse = "")
      
      ### Runs ORA
      system(CommandsToRunORA, wait = T)
    }
  }
}

if (grepl(pattern = ProgramName, ignore.case = T, x = SoftwareToAuc) == T) {
  cat("\n*** Gathering ", ProgramName,  " results ***\n\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      cat("Gathering ", ProgramName, " results: from ", percent, "% of data - iteration ", perm, "\n")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, SuffixPredictions, sep = "", collapse = "")

      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/", ProgramName, "_PERMUTATIONS/", OutdirFromScript, "/", PrefixOutfileLocal, sep="", collapse = "")
      ExpectedOutfileKey     <-paste(ProgramName, percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        DatasetCode<-paste(ProgramName, "__", percent, "__", perm, sep = "", collapse = "")
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, DatasetCode, ListColoursViolin[[ProgramName]], sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Get AUC's from tested software and generate violin plots
####################################

if (grepl(pattern = "Y", ignore.case = T, x = PlotViolin) == T) {
  writeLines("\n*** Get AUC from tested software ***\n")
  
  ListSoftwareToPlot<-unlist(strsplit(SoftwareToAuc, ","))
  
  CommandToGetPerformancePlotsAndAucs<-paste(ListExternalScripts["PathNameToPerformancePlotsAndAucs"],
                                             " -path_outfiles ", Outdir, "/BENCHMARK/MERGED_TABLES/ ",
                                             " -infile_list_infiles ", OutfileListToGetAucs,
                                             " -infile_gold_standards ", InfileGold,
                                             " -use_graphics_device ", "NA",
                                             " -cex_labels ", "NA",
                                             " -lwd  ", "NA",
                                             " -print_plot_ticks_labels_legend  ", "NA",
                                             sep = "", collapse = "")
  system(CommandToGetPerformancePlotsAndAucs, wait = T)
  #
  ExpectedAucROCOutfile  <-paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".PERMUTE.listToGetAuc.ROCcurves.tsv", sep="", collapse = "")
  ExpectedAucPROutfile   <-paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".PERMUTE.listToGetAuc.PRcurves.tsv",  sep="", collapse = "")

  auc.roc.mat<-read.table(file = ExpectedAucROCOutfile, header = T, row.names = NULL)
  auc.pr.mat <-read.table(file = ExpectedAucPROutfile,  header = T, row.names = NULL)
  
  ####################################
  ### Define labels and ylim for violin plots
  ####################################
  writeLines("\n*** Define colours, labels and ylim for violin plots ***\n")
  
  RocAucAxesLimits<-unlist(strsplit(RocAucAxes, ","))
  RocAucMinYLimit <-as.numeric(RocAucAxesLimits[1])
  RocAucMaxYLimit <-as.numeric(RocAucAxesLimits[2])
  RocAucMidY      <-RocAucMinYLimit + ((RocAucMaxYLimit - RocAucMinYLimit) / 2)
  
  PrAucAxesLimits<-unlist(strsplit(PrAucAxes, ","))
  PrAucMinYLimit <-as.numeric(PrAucAxesLimits[1])
  PrAucMaxYLimit <-as.numeric(PrAucAxesLimits[2])
  PrAucMidY      <-PrAucMinYLimit  + ((PrAucMaxYLimit - PrAucMinYLimit) / 2)
  
  YaxisLabelsRoc  <-list()
  YaxisLabelsPR   <-list()
  YaxisLabelsRoc[[1]]<-RocAucMinYLimit
  YaxisLabelsRoc[[2]]<-RocAucMidY
  YaxisLabelsRoc[[3]]<-RocAucMaxYLimit
  YaxisLabelsPR[[1]] <-PrAucMinYLimit
  YaxisLabelsPR[[2]] <-PrAucMidY
  YaxisLabelsPR[[3]] <-PrAucMaxYLimit

  ####################################
  ### Generate violin plots
  ####################################
  writeLines("\n*** Generate violin plots ***\n")
  
  graphics.off()
  plotHeight<-2+length(ListSoftwareToPlot) * PlotHeight
  plotWidth <-length(seq(from = MinPercent, to = MaxPercent, by=PercentBy))
  
  ## It can be called with dev.set(2)
  pdf(paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".ROBUSTNESS_ROC.pdf", sep = "", collapse = ""), width = plotWidth, height = plotHeight)
  par(mfrow=c(length(ListSoftwareToPlot),1), mar= c(3,4,0.1,2) + 0.1)

  ## It can be called with dev.set(3)
  pdf(paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".ROBUSTNESS_PR.pdf",  sep = "", collapse = ""), width = plotWidth, height = plotHeight)
  par(mfrow=c(length(ListSoftwareToPlot),1), mar= c(3,4,0.1,2) + 0.1)
  
  NumberOfProgramPlotted<-0
  for (program in ListSoftwareToPlot) {
    NumberOfProgramPlotted<-NumberOfProgramPlotted+1
    PlotPosition  <-0
    DataToPlotRoc <-list()
    DataToPlotPR  <-list()
    XaxisLabels   <-list()
    ColourForPlots<-paste("#", ListColoursViolin[[program]], sep = "", collapse = "")
    
    for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
      PlotPosition<-PlotPosition+1
      KeyToLookFor<-paste(program, "__", percent, "__", sep = "", collapse = "")
      #
      matchesRoc<-grepl(pattern = KeyToLookFor, x = auc.roc.mat[,"Dataset"])
      DataToPlotRoc[[PlotPosition]] <- as.vector(auc.roc.mat[matchesRoc,"ROC_AUC"])
      #
      matchesPR<-grepl(pattern = KeyToLookFor, x = auc.pr.mat[,"Dataset"])
      DataToPlotPR[[PlotPosition]] <- as.vector(auc.pr.mat[matchesPR,"PR_AUC"])
      #
      XaxisLabels[[PlotPosition]]<-percent
    }

    ### ROC AUC's
    dev.set(2)
    plot(0, type='n', xlim=c(0.5, length(DataToPlotRoc)+0.5), ylim=c(RocAucMinYLimit,RocAucMaxYLimit), xaxt='n', yaxt='n', ylab = "")
    lapply(
      seq_along(DataToPlotRoc), function(percent) {
        vals <- DataToPlotRoc[[percent]]
        if(length(unique(vals)) == 1){
          points(x=percent, y=vals[[1]], pch=16, col = ColourForPlots)
        }else{
          vioplot(DataToPlotRoc[[percent]], at=percent, add=TRUE, col = ColourForPlots)
        }
      }
    )
    if (grepl(pattern = "all|tick_values|tick_val_lab", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      axis(side = 1, at=c(1:length(DataToPlotRoc)), labels = XaxisLabels)
      axis(side = 2, at=c(RocAucMinYLimit,RocAucMidY,RocAucMaxYLimit), labels = YaxisLabelsRoc)
    }else if (grepl(pattern = "ticks", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      axis(side = 1, at=c(1:length(DataToPlotRoc)), labels = FALSE)
      axis(side = 2, at=c(RocAucMinYLimit,RocAucMidY,RocAucMaxYLimit), labels = FALSE)
    }

    if (grepl(pattern = "all|tick_val_lab", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      mtext(text = "ROC  AUC", side=2, line = 2.2)
      if (NumberOfProgramPlotted == length(ListSoftwareToPlot)) {
        mtext(text = "Percentage of genes maintained from original signatures", side=1, line = 2)
      }
    }

    if (grepl(pattern = "all|software", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      text(labels = program,  x = length(XaxisLabels)-1, y = RocAucMinYLimit+0.1)
    }

    abline(h=RocAucAbbline, lty=2, col="gray60", lwd=1)

    ### PR AUC's
    dev.set(3)
    plot(0, type='n', xlim=c(0.5, length(DataToPlotPR)+0.5), ylim=c(PrAucMinYLimit,PrAucMaxYLimit), xaxt='n', yaxt='n', ylab = "")
    lapply(
      seq_along(DataToPlotPR), function(percent) {
        vals <- DataToPlotPR[[percent]]
        if(length(unique(vals)) == 1){
          points(x=percent, y=vals[[1]], pch=16, col = ColourForPlots)
        }else{
          vioplot(DataToPlotPR[[percent]], at=percent, add=TRUE, col = ColourForPlots)
        }
      }
    )

    if (grepl(pattern = "all|tick_values|tick_val_lab", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      axis(side = 1, at=c(1:length(DataToPlotPR)), labels = XaxisLabels)
      axis(side = 2, at=c(PrAucMinYLimit,PrAucMidY,PrAucMaxYLimit), labels = YaxisLabelsPR)
    }else if (grepl(pattern = "ticks", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      axis(side = 1, at=c(1:length(DataToPlotPR)), labels = FALSE)
      axis(side = 2, at=c(PrAucMinYLimit,PrAucMidY,PrAucMaxYLimit), labels = FALSE)
    }

    if (grepl(pattern = "all|tick_val_lab", ignore.case = T, x = PrintLabelsAndTicks) == T) {
      mtext(text = "PR  AUC", side=2, line = 2.2)
      if (NumberOfProgramPlotted == length(ListSoftwareToPlot)) {
        mtext(text = "Percentage of genes maintained from original signatures", side=1, line = 2)
      }
    }

    if (grepl(pattern = "all|software", ignore.case = T, x = PrintLabelsAndTicks) == T) {
    text(labels = program,  x = length(XaxisLabels)-1, y = PrAucMinYLimit+0.1)
    }

    abline(h=PrAucAbbline, lty=2, col="gray60", lwd=1)
  }
  graphics.off()
}

####################################
### Delete temporary files
####################################
writeLines("\n*** Delete temporary files ***\n")
cat(OutfileListToGetAucs)
#file.remove(OutfileListToGetAucs)

####################################
### Report time used
####################################
writeLines("\n*** Report time used ***\n")

EndTimeOverall<-Sys.time()

TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "secs"))

OutfileCPUusage<-paste(Outdir,"/",PrefixOutfiles,".PERMUTE_CPUusage.txt", sep="")
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

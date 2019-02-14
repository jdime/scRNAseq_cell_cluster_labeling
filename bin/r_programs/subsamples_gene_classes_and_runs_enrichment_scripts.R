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
### (i)   Formats inputs for each requested cell type prediction method: CIBERSORT, GSEA, GSVA, and ORA
### (ii)  Subsamples cell type gene expression signatures using infiles (2)
### (iii) Runs cell type predictions methods using infile (1) and subsampled signatures from step (ii)
### (iv)  Gather results from step (iii)
### (v)   Runs ROC, PR and AUC analyses using results from step (iv) and infile (3)
### (vi)  Runs robustness analyses (optional), including violin plots of AUC distributions with results from step (v)
###
### Notes:
###        Note 1
###        In the case of CIBERSORT, two forms of cell type signature infiles can be provided:
###        - in the form of gene sets (2a) or in the form of gene expression profiles (2b)
###        - Subsampling genes from is always conducted using 2a type files, which can be propagated to 2b type files using script:
###          'propagates_permuted_gmt_files_to_profile.R'
###
###        Note 2
###        To add new methods to evaluate you can:
###        a) add the path to the executable in 'Dependencies:'
###        b) provide inputs and parameters adding commands indicated by 'Runs NEW_METHOD and gather its outputs'
###        c) add program name to this script 'option_list' parameters --software_to_run and --software_to_auc
###        d) add program name to script 'obtains_performance_plots_from_cluster_labelings.pl' so that it can take a 'Generic_ES' infile
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
### 'optparse'   (CRAN) to handle one-line-commands
### 'vioplot'    (CRAN) to generate violin plots of permutation results
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(vioplot))
###
### External scripts (check each script for their own dependencies)
### and change their path as needed here:
PathNameToPermuteGmtRscript       <- "~/r_programs/obtains_permuted_samples_from_gmt.R"
PathNameToPermuteProfileRscript   <- "~/r_programs/propagates_permuted_gmt_files_to_profile.R"
PathNameToCibersortscript         <- "~/bin/obtains_CIBERSORT_for_MatrixColumns.pl"
PathNameToGSEAscript              <- "~/bin/obtains_GSEA_for_MatrixColumns.pl"
PathNameToGSVAscript              <- "~/r_programs/obtains_GSVA_for_MatrixColumns.R"
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
              ...etc"),
  
  make_option(c("-c", "--infile_signature_gmt"), default="NA",
              help="Path/name to a <tab> delimited cell type signature *file* in the form of gene lists (*gmt format), like:
              GeneSet1_ID  GeneSet1_Name  Gene1 Gene2 ... etc
              GeneSet2_ID  GeneSet2_Name  Gene2 Gene3 ... etc"),
  
  make_option(c("-g", "--infile_gold"), default="NA",
              help="A path/name to a <tab> delimited *file* of gold standard cluster labels in format, like:
              clust1  GeneSet1_ID
              clust2  GeneSet3_ID
              clust3  GeneSet5_ID
              ... etc"),
  
  make_option(c("-t", "--permute_gmt"), default="n",
              help="Indicates if permutation of --infile_signature_gmt should be conducted [use y/Y] or permuted files exist already [type n/N]. Default 'n'."),
  
  make_option(c("-u", "--propagate_permuted_gmt_to_profiles"), default="n",
              help="Only needed by 'CIBERSORT_PROFILE'
              Indicates if permutation of --infile_signature_gmt files should be propagated to profiles [use y/Y] or permuted files exist already [type n/N]. Default 'n'."),

  make_option(c("-d", "--infile_signature_profile"), default="NA",
              help="Only needed by 'CIBERSORT_PROFILE'
              Path/name to a <tab> delimited cell type signature *file* in the form of gene expression profile matrix, like:
              GENE   CellType1  CellType2  CellType3 ... etc
              Gene1  0.1        0.04       0.6
              Gene2  0.02       0.1        0.01
              Gene3  0.04       0.3        0.06
              ...etc"),

  make_option(c("-s", "--sample_percentages"), default="10,100,10",
              help="Indicates three values: minimum_percentage, maximum_percentage, increment to be used to sample --infile_signature_gmt. Default '10,100,10'.
              For example, to subsample percentages from 10 up to 50, increasing by 5 (i. e. 10, 15, 20, ... 50) use '10,50,5'
              Note: if using '-n N' then this script will take previously permuted *gmt infiles, stored at /path_to/--outdir/GMT_PERMUTATIONS"),
  
  make_option(c("-r", "--iterations"), default="1-100",
              help="Indicates the number of times to subsample --infile_signature_gmt, e.g. '1-100' to run from iteration 1 to 100,
              or '50-100' to run from iteration 50 to 100. Default '1-100'.
              Note: if using '-t N', then this script will use previously permuted *gmt infiles"),
  
  make_option(c("-m", "--software_to_run"), default="NONE",
              help="Indicates <comma> delimited name(s) of enrichment sofware to *run*:
              CIBERSORT_BINARY,CIBERSORT_PROFILE,GSEA,GSVA,ORA
              Default 'NONE'."),
  
  make_option(c("-n", "--software_to_auc"), default="NONE",
              help="Indicates <comma> delimited name(s) of enrichment sofware results to *gather* together:
              CIBERSORT_BINARY,CIBERSORT_PROFILE,GSEA,GSVA,ORA
              This may be useful to gather results from separate runs
              Default 'NONE'."),
  
  make_option(c("-v", "--generate_violin_plots"), default="n",
              help="Indicates if violin plots should be generated from --software_to_run and --software_to_auc steps. Type [y/Y] or [n/N]. Default 'n'."),
  
  make_option(c("-a", "--roc_auc_violin_plot_y_axes_limits"), default="0.5,1",
              help="Indicates the min,max values for the ROC AUC violin plot y-axes of enrichment software. Default '0.5,1'."),
  
  make_option(c("-b", "--pr_auc_violin_plot_y_axes_limits"), default="0,1",
              help="Indicates the min,max values for the PR AUC violin plot y-axes of enrichment software. Default '0,1'."),
  
  make_option(c("-k", "--roc_horizontal_guide_line"), default="NA",
              help="Indicates a value in the y-axis of the ROC AUC violin plots to place a horizontal guide line. Or type 'NA' to omit. Default 'NA'."),
  
  make_option(c("-l", "--pr_horizontal_guide_line"), default="NA",
              help="Indicates a value in the y-axis of the PR AUC violin plots to place a horizontal guide line. Or type 'NA' to omit. Default 'NA'."),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
              Note this script will automatically add 'percentXX.permYY' to the outfile name indicating -s value and iteration number. It can't be 'NA'."),
  
  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory. It can't be 'NA'."),
  
  make_option(c("-y", "--ora_use_values_or_ranks"), default="NA",
              help="This applies only to the ORA test:
              Indicates if a cutoff applied to --infile_mat to get 'white/drawn' balls should be based on scores 'ranks' or 'values'
              For example, using '-y values -z 1.5' will use all genes with values >= 1.5 in each column of --infile_mat as 'drawn' balls
              Whereas, using '-y ranks -z 100' will sort each column of --infile_mat and use the top 100 ranked genes as 'drawn' balls
              Or type 'NA' if ORA method is not being used"),
  
  make_option(c("-z", "--ora_mat_cutoff"), default="NA",
              help="This applies only to the ORA test:
              Indicates the cutoff to be applied to --infile_mat to get 'white/drawn' balls
              Or type 'NA' if ORA method is not being used")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat         <- opt$infile_mat
InfileGmt         <- opt$infile_signature_gmt
InfileProfile     <- opt$infile_signature_profile
InfileGold        <- opt$infile_gold
Outdir            <- opt$outdir
PrefixOutfiles    <- opt$prefix_outfiles
Iterations        <- opt$iterations
SamplePercentages <- opt$sample_percentages
OraValOrRank      <- opt$ora_use_values_or_ranks
OraCutoff         <- opt$ora_mat_cutoff
SoftwareToRun     <- opt$software_to_run
SoftwareToAuc     <- opt$software_to_auc
PermuteGmt        <- opt$permute_gmt
PermuteProfile    <- opt$propagate_permuted_gmt_to_profiles
PlotViolin        <- opt$generate_violin_plots
RocAucAxes        <- opt$roc_auc_violin_plot_y_axes_limits
PrAucAxes         <- opt$pr_auc_violin_plot_y_axes_limits
RocAucAbbline     <- opt$roc_horizontal_guide_line
PrAucAbbline      <- opt$pr_horizontal_guide_line

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_mat", "infile_signature_gmt", "infile_gold", "outdir", "prefix_outfiles", "ora_use_values_or_ranks", "ora_mat_cutoff")
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
  print("Not permutting --infile_signature_gmt. Assuming permuted *gmt files already exist")
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
    stop(paste("Couldn't find file: '", InfileGmtLocal,  "'",sep="", collapse = ""))
  }
}else if (grepl(pattern = "CIBERSORT_PROFILE", ignore.case = T, x = SoftwareToRun) == T) {
  print("Not permutting --infile_signature_profile. Assuming permuted gene expression profile files already exist")
}

####################################
### Runs CIBERSORT_BINARY and gather its outputs
####################################
if (grepl(pattern = "CIBERSORT_BINARY", ignore.case = T, x = SoftwareToRun) == T) {
  writeLines("\n*** Runing CIBERSORT_BINARY ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Running CIBERSORT_BINARY: with ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      #
      InfileGmtLocal<-paste(Outdir,"/GMT_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunCIBERSORT_BINARY<-paste(ListExternalScripts["PathNameToCibersortscript"],
                                           " -infile_matrix ",   InfileMat,
                                           " -infile_classes ",  InfileGmtLocal,
                                           " -path_outfiles ",   Outdir, "/CIBERSORT_BINARY_PERMUTATIONS/",
                                           " -prefix_outfiles ", PrefixOutfileLocal,
                                           " -nperm ",           1000,
                                           " -classes_type gmt ",
                                           sep = "", collapse = "")
      
      ### Runs CIBERSORT_BINARY
      system(CommandsToRunCIBERSORT_BINARY, wait = T)
    }
  }
}

if (grepl(pattern = "CIBERSORT_BINARY", ignore.case = T, x = SoftwareToAuc) == T) {
  writeLines("\n*** Gathering CIBERSORT_BINARY results ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Gathering CIBERSORT_BINARY results: from ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/CIBERSORT_BINARY_PERMUTATIONS/CIBERSORT/", PrefixOutfileLocal, ".CIBERSORT_enrichment_scores.tsv", sep="", collapse = "")
      ExpectedOutfileKey     <-paste("CIBERSORT_BINARY", percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, paste("CIBERSORT_BINARY", "__", percent, "__", perm, sep = "", collapse = ""), "Generic_ES", sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs CIBERSORT_PROFILE and gather its outputs
####################################
if (grepl(pattern = "CIBERSORT_PROFILE", ignore.case = T, x = SoftwareToRun) == T) {
  writeLines("\n*** Runing CIBERSORT_PROFILE ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Running CIBERSORT_PROFILE: with ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      #
      InfileProfileLocal<-paste(Outdir,"/PROFILE_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".profile.tsv", sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunCIBERSORT_PROFILE<-paste(ListExternalScripts["PathNameToCibersortscript"],
                                            " -infile_matrix ",   InfileMat,
                                            " -infile_classes ",  InfileProfileLocal,
                                            " -path_outfiles ",   Outdir, "/CIBERSORT_PROFILE_PERMUTATIONS/",
                                            " -prefix_outfiles ", PrefixOutfileLocal,
                                            " -nperm ",           1000,
                                            " -classes_type profile ",
                                            sep = "", collapse = "")
      
      ### Runs CIBERSORT_PROFILE
      system(CommandsToRunCIBERSORT_PROFILE, wait = T)
    }
  }
}

if (grepl(pattern = "CIBERSORT_PROFILE", ignore.case = T, x = SoftwareToAuc) == T) {
  writeLines("\n*** Gathering CIBERSORT_PROFILE results ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Gathering CIBERSORT_PROFILE results: from ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/CIBERSORT_PROFILE_PERMUTATIONS/CIBERSORT/", PrefixOutfileLocal, ".CIBERSORT_enrichment_scores.tsv", sep="", collapse = "")
      ExpectedOutfileKey     <-paste("CIBERSORT_PROFILE", percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, paste("CIBERSORT_PROFILE", "__", percent, "__", perm, sep = "", collapse = ""), "Generic_ES", sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs GSEA and gather its outputs
####################################
if (grepl(pattern = "GSEA", ignore.case = T, x = SoftwareToRun) == T) {
  writeLines("\n*** Runing GSEA ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Running GSEA: with ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      #
      InfileGmtLocal<-paste(Outdir,"/GMT_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunGSEA<-paste(ListExternalScripts["PathNameToGSEAscript"],
                               " -infile_matrix ",   InfileMat,
                               " -infile_classes ",  InfileGmtLocal,
                               " -path_outfiles ",   Outdir, "/GSEA_PERMUTATIONS/",
                               " -prefix_outfiles ", PrefixOutfileLocal,
                               " -cutoff_print_p ",  0.01,
                               " -cutoff_print_q ",  0.25,
                               " -nperm ",           1000,
                               sep = "", collapse = "")
      
      ### Runs GSEA
      system(CommandsToRunGSEA, wait = T)
    }
  }
}

if (grepl(pattern = "GSEA", ignore.case = T, x = SoftwareToAuc) == T) {
  writeLines("\n*** Gathering GSEA results ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Gathering GSEA results: from ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir, "/GSEA_PERMUTATIONS/GSEA/", PrefixOutfileLocal, "/ClassLevel/", PrefixOutfileLocal, ".Pval.Unfiltered.mat.txt", sep="", collapse = "")
      ExpectedOutfileKey     <-paste("GSEA", percent, perm, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, paste("GSEA", "__", percent, "__", perm, sep = "", collapse = ""), "GSEA_Pvalues", sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs GSVA and gather its outputs
####################################
if (grepl(pattern = "GSVA", ignore.case = T, x = SoftwareToRun) == T) {
  writeLines("\n*** Runing GSVA ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Running GSVA: with ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      #
      InfileGmtLocal<-paste(Outdir,"/GMT_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      CommandsToRunGSVA<-paste("Rscript ", ListExternalScripts["PathNameToGSVAscript"],
                               " -i ", InfileMat,
                               " -c ", InfileGmtLocal,
                               " -o ", Outdir, "/GSVA_PERMUTATIONS/",
                               " -p ", PrefixOutfileLocal,
                               " -e 0.05 -f 0.1",
                               sep = "", collapse = "")
      ### Runs GSVA
      system(CommandsToRunGSVA, wait = T)
    }
  }
}

if (grepl(pattern = "GSVA", ignore.case = T, x = SoftwareToAuc) == T) {
  writeLines("\n*** Gathering GSVA results ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Gathering GSVA results: from ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/GSVA_PERMUTATIONS/GSVA/", PrefixOutfileLocal, ".GSVA_all_scores_table.tsv", sep="", collapse = "")
      ExpectedOutfileKey     <-paste("GSVA", percent, perm, sep="\t", collapse = "\t")
      
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, paste("GSVA", "__", percent, "__", perm, sep = "", collapse = ""), "GSVA_ES", sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Runs ORA and gather its outputs
####################################
if (grepl(pattern = "ORA", ignore.case = T, x = SoftwareToRun) == T) {
  writeLines("\n*** Runing ORA ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Running ORA: with ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      #
      InfileGmtLocal<-paste(Outdir,"/GMT_PERMUTATIONS/", PrefixOutfiles, ".percent", percent, ".perm", perm, ".gmt", sep="", collapse = "")
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, ".", OraValOrRank, OraCutoff, sep = "", collapse = "")
      CommandsToRunORA<-paste(ListExternalScripts["PathNameToORAscript"],
                                 " -infile_matrix ",       InfileMat,
                                 " -infile_classes ",      InfileGmtLocal,
                                 " -path_outfiles ",       Outdir, "/ORA_PERMUTATIONS/",
                                 " -prefix_outfiles ",     PrefixOutfileLocal,
                                 " -cutoff_pos ",          OraCutoff,
                                 " -use_values_or_ranks ", OraValOrRank,
                                 " -cutoff_ora ",          "0.05",
                                 " -cutoff_neg  ",         "NA",
                                 " -numb_cols ",           "ALL",
                                 " -p_correction_test ",   "BH",
                                 " -use_ora ",             "pc",
                                 " -restrict_classes ",    "ALL",
                                 " -set_max ",             "1000",
                                 " -set_min ",             "1",
                                 " -use_universe ",        "i",
                                 sep = "", collapse = "")
      
      ### Runs ORA
      system(CommandsToRunORA, wait = T)
    }
  }
}

if (grepl(pattern = "ORA", ignore.case = T, x = SoftwareToAuc) == T) {
  writeLines("\n*** Gathering ORA results ***\n")
  
  for (percent in seq(from = MinPercent, to = MaxPercent, by=PercentBy)) {
    for (perm in c(MinIterations:MaxIterations)) {
      print(paste("Gathering ORA results: from ", percent, "% of data - iteration ", perm, sep = "", collapse = ""))
      PrefixOutfileLocal<-paste(PrefixOutfiles, ".percent", percent, ".perm", perm, ".", OraValOrRank, OraCutoff, sep = "", collapse = "")
      
      ### Adds expected outfile to list for compilation
      ExpectedOutfilePathName<-paste(Outdir,"/ORA_PERMUTATIONS/ORA/", PrefixOutfileLocal, ".ORA.PvalUncorrected.All.mat", sep="", collapse = "")
      ExpectedOutfileKey     <-paste("ORA", percent, perm, ".", OraValOrRank, OraCutoff, sep="\t", collapse = "\t")
      if(file.exists(ExpectedOutfilePathName)) {
        ListExpectedOutfiles[ExpectedOutfileKey]<-ExpectedOutfilePathName
        write(file = OutfileListToGetAucs, append = T, x=paste(ExpectedOutfilePathName, paste("ORA", "__", percent, "__", perm, sep = "", collapse = ""), "ORA_Pvalues", sep = "\t", collapse = "\t"))
      }else{
        stop(paste("Couldn't find file: '", ExpectedOutfilePathName,  "'",sep="", collapse = ""))
      }
    }
  }
}

####################################
### Get AUC's from tested software
####################################
if (grepl(pattern = "Y", ignore.case = T, x = PlotViolin) == T) {
  writeLines("\n*** Get AUC from tested software ***\n")
  
  ListSoftwareToPlot<-unlist(strsplit(SoftwareToAuc, ","))
  
  CommandToGetPerformancePlotsAndAucs<-paste(ListExternalScripts["PathNameToPerformancePlotsAndAucs"],
                                             " -path_outfiles ", Outdir, "/BENCHMARK/MERGED_TABLES/ ",
                                             " -infile_list_infiles ", OutfileListToGetAucs,
                                             " -infile_gold_standards ", InfileGold,
                                             " -use_graphics_device ", "NA",
                                             sep = "", collapse = "")
  system(CommandToGetPerformancePlotsAndAucs, wait = T)
  #
  ExpectedAucROCOutfile  <-paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".PERMUTE.listToGetAuc.ROCcurves.tsv", sep="", collapse = "")
  ExpectedAucPROutfile   <-paste(Outdir, "/BENCHMARK/MERGED_TABLES/PERFORMANCE_PLOTS/", PrefixOutfiles, ".PERMUTE.listToGetAuc.PRcurves.tsv",  sep="", collapse = "")
  auc.roc.mat<-read.table(file = ExpectedAucROCOutfile, header = T, row.names = NULL)
  auc.pr.mat <-read.table(file = ExpectedAucPROutfile,  header = T, row.names = NULL)
  
  ####################################
  ### Define colours, labels and ylim for violin plots
  ####################################
  writeLines("\n*** Define colours, labels and ylim for violin plots ***\n")
  
  ListColoursViolin<-list("CIBERSORT_BINARY"  = "#E69F00", # reddishpurple
                          "CIBERSORT_PROFILE" = "#009E73", # vermillion
                          "GSEA"              = "#CC79A7", # bluishgreen
                          "GSVA"              = "#1E90FF", # skyblue
                          "ORA"            = "#D55E00"  # orange
  )
  
  RocAucAxesLimits<-unlist(strsplit(RocAucAxes, ","))
  RocAucMinYLimit <-as.numeric(RocAucAxesLimits[1])
  RocAucMaxYLimit <-as.numeric(RocAucAxesLimits[2])

  PrAucAxesLimits<-unlist(strsplit(PrAucAxes, ","))
  PrAucMinYLimit <-as.numeric(PrAucAxesLimits[1])
  PrAucMaxYLimit <-as.numeric(PrAucAxesLimits[2])
  
  ####################################
  ### Generate violin plots
  ####################################
  writeLines("\n*** Generate violin plots ***\n")
  
  graphics.off()
  
  pdf(paste(Outdir,"/",PrefixOutfiles,".ROBUSTNESS_ROC.pdf", sep = "", collapse = ""), width = 7, height = 8) ## It can be called with dev.set(2)
  par(mfrow=c(length(ListSoftwareToPlot),1), mar= c(3,4,0.1,2) + 0.1)

  pdf(paste(Outdir,"/",PrefixOutfiles,".ROBUSTNESS_PR.pdf",  sep = "", collapse = ""), width = 7, height = 8) ## It can be called with dev.set(2)
  par(mfrow=c(length(ListSoftwareToPlot),1), mar= c(3,4,0.1,2) + 0.1)
  
  for (program in ListSoftwareToPlot) {
    PlotPosition<-0
    DataToPlotRoc<-list()
    DataToPlotPR <-list()
    XaxisLabels  <-list()
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
    plot(0, type='n', xlim=c(0.5, length(DataToPlotRoc)+0.5), ylim=c(RocAucMinYLimit,RocAucMaxYLimit), xaxt='n', ylab = "")
    lapply(seq_along(DataToPlotRoc), function(percent) vioplot(DataToPlotRoc[[percent]], at=percent, add=TRUE, col = ListColoursViolin[[program]]))
    axis(side = 1, at=c(1:length(DataToPlotRoc)), labels = XaxisLabels)
    text(labels = program,  x = length(XaxisLabels)-1, y = RocAucMinYLimit+0.1)
    abline(h=RocAucAbbline, lty=2, col="gray60", lwd=1)
    mtext(text = "ROC  AUC", side=2, line = 2.2, cex=0.75)

    ### PR AUC's
    dev.set(3)
    plot(0, type='n', xlim=c(0.5, length(DataToPlotPR)+0.5), ylim=c(PrAucMinYLimit,PrAucMaxYLimit), xaxt='n', ylab = "")
    lapply(seq_along(DataToPlotPR), function(percent) vioplot(DataToPlotPR[[percent]], at=percent, add=TRUE, col = ListColoursViolin[[program]]))
    axis(side = 1, at=c(1:length(DataToPlotPR)), labels = XaxisLabels)
    text(labels = program,  x = length(XaxisLabels)-1, y = PrAucMinYLimit+0.1)
    abline(h=PrAucAbbline, lty=2, col="gray60", lwd=1)
    mtext(text = "PR  AUC", side=2, line = 2.2, cex=0.75)
    
  }
  graphics.off()
}

####################################
### Delete temporary files
####################################
writeLines("\n*** Delete temporary files ***\n")

file.remove(OutfileListToGetAucs)

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
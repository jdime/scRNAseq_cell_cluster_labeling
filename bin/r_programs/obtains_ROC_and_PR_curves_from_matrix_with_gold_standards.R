####################################
### Script 'obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R' obtains Receiver Operating Characteristic (ROC)
### and Precision-Recall (PR) curves for each column of predictions in a --infile_mat
### --infile_mat must have binary gold standard labels in column #2, and predictions in columns 3 to n (see 'option_list' below)
###
### Notes: two types of PR curves are obtained 'saw-toothed' and 'interpolated'.
###        Differences between the two types are documented here:
###        https://rdrr.io/cran/DMwR/man/PRcurve.html  and
###        https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-ranked-retrieval-results-1.html
####################################
### Questions/comments to Javier Diaz - javier.diazmejia@gmail.com
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to_this_file/obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R -h'
### for help
####################################

####################################
### Dependencies:
####################################
suppressPackageStartupMessages(library(optparse)) # (CRAN) to handle one-line-commands
suppressPackageStartupMessages(library(precrec))  # (CRAN) to generate ROC and 'saw-toothed' PR curves, and to get ROC AUC and PR AUC
suppressPackageStartupMessages(library(ROCR))     # (CRAN) to generate interpolated PR curves
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
              help="A path/name to a <tab> delimited *file* with two or more columns with predictions values
              and another with binary labels ('1' for positive, '0' for negative gold standards):
              ID_row   Label  Prediction1  Prediction2
              row_1    1      0.612        0.612
              row_2    0      0.364        0.364
              row_3    1      0.432        0.432
              row_4    0      0.140        0.140"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID"),

  make_option(c("-l", "--lwd"), default="3",
              help="The line width, a positive number, defaulting to 3"),

  make_option(c("-d", "--graphics_device"), default="pdf",
              help="Indicates if the plots should be either 'pdf' (default) or 'png'
              or type 'NA' to only generate AUC *.tsv files")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
Lwd             <- opt$lwd
GraphicsDevice  <- opt$graphics_device
Tempdir         <- "~/temp" ## Using this for temporary storage of outfiles because sometimes long paths of outdirectories casuse R to leave outfiles unfinished

StartTimeOverall<-Sys.time()

####################################
### Check that mandatory parameters are not 'NA' (default)
####################################

ListMandatory<-list("infile_mat")
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
Outdir<-gsub("^~/", paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Outdir)
Tempdir<-gsub("^~/",paste(c(UserHomeDirectory,"/"), sep = "", collapse = ""), Tempdir)
Outdir<-gsub("/$", "", Outdir)
Tempdir<-gsub("/$", "", Tempdir)
#
dir.create(file.path(Outdir),  showWarnings = F, recursive = T)
dir.create(file.path(Tempdir), showWarnings = F, recursive = T)

OutRocTsv  <-paste(Tempdir, "/", PrefixOutfiles, ".ROCcurves.tsv", sep = "", collapse = "")
OutPrStTsv <-paste(Tempdir, "/", PrefixOutfiles, ".PRcurves.tsv", sep = "", collapse = "")

graphics.off()
GeneratePlots<-0
if (regexpr("^pdf$", GraphicsDevice, ignore.case = T)[1] == 1) {
  GeneratePlots<-1
  OutRocPlot  <-paste(Tempdir, "/", PrefixOutfiles, ".ROCcurves.pdf", sep = "", collapse = "")
  OutPrStPlot <-paste(Tempdir, "/", PrefixOutfiles, ".PR_saw_tooth_curves.pdf", sep = "", collapse = "")
  OutPrItPlot <-paste(Tempdir, "/", PrefixOutfiles, ".PR_interpolated_curves.pdf", sep = "", collapse = "")
  pdf(file = OutRocPlot,  width = 7, height = 7) ## It can be called with dev.set(2)
  pdf(file = OutPrStPlot, width = 7, height = 7) ## It can be called with dev.set(3)
  pdf(file = OutPrItPlot, width = 7, height = 7) ## It can be called with dev.set(4)
}else if (regexpr("^png$", GraphicsDevice, ignore.case = T)[1] == 1) {
  GeneratePlots<-1
  OutRocPlot <-paste(Tempdir, "/", PrefixOutfiles, ".ROCcurves.png", sep = "", collapse = "")
  OutPrStPlot  <-paste(Tempdir, "/", PrefixOutfiles, ".PR_saw_tooth_curves.png", sep = "", collapse = "")
  OutPrItPlot <-paste(Tempdir, "/", PrefixOutfiles, ".PR_interpolated_curves.png", sep = "", collapse = "")
  png(file = OutRocPlot,  width = 480, height = 480) ## It can be called with dev.set(2)
  png(file = OutPrStPlot, width = 480, height = 480) ## It can be called with dev.set(3)
  png(file = OutPrItPlot, width = 480, height = 480) ## It can be called with dev.set(4)
}else if (regexpr("^NA$", GraphicsDevice, ignore.case = T)[1] == 1) {
  print("Only text *tsv files will be generated with AUC's")
}else{
  stop(paste("Unexpected type of graphics device: ", GraphicsDevice, "\n\nFor help type:\n\nRscript obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R -h\n\n", sep=""))
}

####################################
### Define colours
####################################
writeLines("\n*** Define colours ***\n")

ColoursNumbToHex <- list(
  "1"  = "#E69F00",
  "2"  = "#009E73",
  "3"  = "#CC79A7",
  "4"  = "#1E90FF",
  "5"  = "#D55E00",
  "6"  = "#8B8989",
  "7"  = "#FFD700",
  "8"  = "#43CD80",
  "9"  = "#FFC0CB",
  "10" = "#56B4E9",
  "11" = "#8B4789",
  "12" = "#0072B2",
  "13" = "#000000"
)

####################################
### Load data, get plots and AUC's
####################################
writeLines("\n*** Load data, get plots and AUC's ***\n")

mat<-read.table(InfileMat, header = T, row.names = 1, sep = "\t")

####################################
### Generate plots and get AUC's
####################################
writeLines("\n*** Generate plots and get AUC's ***\n")

labelsToPlot<-mat[,1]

AUCs_ROC_vals <-list()
AUCs_PR_vals  <-list()
AUCs_ROC_text <-list()
AUCs_PR_text  <-list()
AUCs_ROC_cols <-list()
AUCs_PR_cols  <-list()
LineTypes     <-list() 

for (i in 2:length(colnames(mat))) {
  colheaderPrediction<-colnames(mat[i])
  print(paste(colnames(mat[1]), " vs. ", colnames(mat[i]), collapse = "", sep = ""))
  predsToPlot<-mat[,i]
  
  ## Get ROC and PR 'saw-toothed' curve coordinates and AUC values
  curves.r <- evalmod(scores = predsToPlot, labels = labelsToPlot)
  aucs.r   <- auc(curves.r)
  AUCs_ROC_vals[i] <-round(aucs.r$aucs[[1]], digits = 2)
  AUCs_PR_vals[i]  <-round(aucs.r$aucs[[2]], digits = 2)
  AUCs_ROC_text[i]<-paste(colnames(mat)[[i]], "=", (round(aucs.r$aucs[[1]], digits = 2)), sep = "", collapse = "")
  AUCs_PR_text[i] <-paste(colnames(mat)[[i]], "=", (round(aucs.r$aucs[[2]], digits = 2)), sep = "", collapse = "")
  AUCs_ROC_cols[i]<-ColoursNumbToHex[i-1][[1]]
  AUCs_PR_cols[i] <-ColoursNumbToHex[i-1][[1]]
  ### This forces lty to start from 2, because lty=1 is a solid line,
  ### which doesn't allow to see overlapping lines
  LineTypes[i] <- as.integer(i / length(ColoursNumbToHex)) + 2
  
  ## Get PR 'interpolated' curve coordinates
  pred <- prediction(predsToPlot,labelsToPlot)
  perf <- performance(pred,"prec","sens")
  x.values<-unlist(slot(perf,"x.values"))[is.finite(unlist(slot(perf,"x.values")))]
  y.values<-unlist(slot(perf,"y.values"))[is.finite(unlist(slot(perf,"y.values")))]
  length(x.values)
  length(y.values)
  y.values.inv<-y.values[length(y.values):1]
  cm<-cummax(y.values.inv[1:length(y.values.inv)])

  if (GeneratePlots == 1) {

    ## Plot ROC curves
    dev.set(2)
    if (i == 2) {
      plot(x=curves.r$rocs[[1]]$x, y=curves.r$rocs[[1]]$y, cex.lab=1.4, cex.axis=1.4, lwd=Lwd, col=AUCs_ROC_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab = "Sensitivity")
      par(new=T)
      plot(c(1:0), c(1:0), lty=2, col="gray60", lwd=1, type="l", axes=FALSE, bty="n", xlab="", ylab="", main="")
    }else{
  		par(new=T)
      plot(x=curves.r$rocs[[1]]$x, y=curves.r$rocs[[1]]$y, cex.lab=1.4, cex.axis=1.4, lwd=Lwd, col=AUCs_ROC_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }

    ## Plot PR 'saw-toothed' curves
    dev.set(3)
    if (i == 2) {
      plot(x=curves.r$prcs[[1]]$x, y=curves.r$prcs[[1]]$y, cex.lab=1.4, cex.axis=1.4, lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xlab = "Recall", ylab = "Precision")
      abline(h=0.5, lty=2, col="gray60", lwd=1)
    }else{
      par(new=T)
      plot(x=curves.r$prcs[[1]]$x, y=curves.r$prcs[[1]]$y, cex.lab=1.4, cex.axis=1.4, lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }
    
    ## Plot PR 'interpolated' curves
    dev.set(4)
    if (i == 2) {
      plot(x=x.values[2:length(x.values)],y=cm[length(cm):1],cex.lab=1.4,cex.axis=1.4,lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xlab = "Recall", ylab = "Precision")
      abline(h=0.5, lty=2, col="gray60", lwd=1)
    }else{
      par(new=T)
      plot(x=x.values[2:length(x.values)],y=cm[length(cm):1], cex.lab=1.4, cex.axis=1.4, lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }
  }
}

if (GeneratePlots == 1) {
  ### Add legends
  ColorsListROC <-as.character(c(AUCs_ROC_cols[2:length(colnames(mat))]))
  ColorsListPR  <-as.character(c(AUCs_ROC_cols[2:length(colnames(mat))]))
  LineTypesList <-as.numeric(c(LineTypes[2:length(colnames(mat))]))

  dev.set(2)
  legend("bottomright", legend=c(AUCs_ROC_text[2:length(colnames(mat))]), col=ColorsListROC, lty=LineTypesList, pch="", cex=1, ncol=1, bty="n",lwd=Lwd)
  dev.set(3)
  legend("bottomleft",  legend=c(AUCs_PR_text[2:length(colnames(mat))]),  col=ColorsListPR,  lty=LineTypesList, pch="", cex=1, ncol=1, bty="n",lwd=Lwd)
  dev.set(4)
  legend("bottomleft",  legend=c(AUCs_PR_text[2:length(colnames(mat))]),  col=ColorsListPR,  lty=LineTypesList, pch="", cex=1, ncol=1, bty="n",lwd=Lwd)
  
  ### Close graphic devices
  graphics.off()
}

####################################
### Generate AUC's outfiles
####################################
writeLines("\n*** Generate AUC's outfiles ***\n")

write(x=paste("Dataset", "ROC_AUC", sep = "\t", collapse = ""), file = OutRocTsv)
write(x=paste(colnames(mat)[2:length(colnames(mat))], AUCs_ROC_vals[2:length(colnames(mat))], "\n", sep = "\t", collapse = ""), append = T, file = OutRocTsv)

write(x=paste("Dataset", "PR_AUC", sep = "\t", collapse = ""), file = OutPrStTsv)
write(x=paste(colnames(mat)[2:length(colnames(mat))], AUCs_PR_vals[2:length(colnames(mat))], "\n", sep = "\t", collapse = ""), append = T, file = OutPrStTsv)

####################################
### Report used options
####################################
writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"/",PrefixOutfiles,".ROC_and_PR_curves_UsedOptions.txt", sep="")
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

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".ROC_and_PR_curves_CPUusage.txt", sep="")
ReportTime<-c(
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

####################################
### Moving outfiles into outdir
####################################
writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir))

outfiles_to_move <- list.files(Tempdir, pattern = paste(PrefixOutfiles, ".ROC_and_PR_curves_", sep=""), full.names = F)
outfiles_to_move <- c(outfiles_to_move, (list.files(Tempdir, pattern = paste(PrefixOutfiles, ".ROCcurves", sep=""), full.names = F)))
outfiles_to_move <- c(outfiles_to_move, (list.files(Tempdir, pattern = paste(PrefixOutfiles, ".PRcurves",  sep=""), full.names = F)))
outfiles_to_move <- c(outfiles_to_move, (list.files(Tempdir, pattern = paste(PrefixOutfiles, ".PR_saw_tooth_curves",  sep=""), full.names = F)))
outfiles_to_move <- c(outfiles_to_move, (list.files(Tempdir, pattern = paste(PrefixOutfiles, ".PR_interpolated_curves",  sep=""), full.names = F)))
outfiles_to_move
sapply(outfiles_to_move,FUN=function(eachFile){
  ### using two steps instead of just 'file.rename' to avoid issues with path to ~/temp in cluster systems
  file.copy(from=paste(Tempdir,"/",eachFile,sep=""),to=paste(Outdir, "/", eachFile, sep=""),overwrite=T)
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

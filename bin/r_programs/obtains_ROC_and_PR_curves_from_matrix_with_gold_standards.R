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
                row_4    0      0.140        0.140
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-c", "--infile_colours"), default="NA",
              help="A path/name to a <tab> delimited *file* with HEX colours for each prediction column from --infile_mat, like:
                Prediction1  #E69F00
                Prediction2  #009E73
                Prediction3  #CC79A7
                ...etc
                Type 'NA' to use up to 13 default colours, defined by list(ColoursNumbToHex)
                Default = 'NA'"),

  make_option(c("-o", "--outdir"), default="NA",
              help="A path/name for the results directory
                Default = 'No default. It's mandatory to specify this parameter'"),
  
  make_option(c("-p", "--prefix_outfiles"), default="NA",
              help="A prefix for outfile names, e.g. your project ID
                Default = 'No default. It's mandatory to specify this parameter'"),

  make_option(c("-l", "--lwd"), default="3",
              help="Indicates the width of plot lines
                Default = '3'"),
  
  make_option(c("-w", "--print_plot_ticks_labels_legend"), default="all",
              help="Indicates if plots should contain:
                a) axes labels, tick marks and legends (type 'all')
                b) only tick marks (type 'ticks')
                c) only axes labels and tick marks (type 'label_tick')
                d) only legend (type 'legend')
                e) no axes labels, tick marks or legends (type 'none')
                Default = 'all'"),
  
  make_option(c("-x", "--cex_labels"), default="1",
              help="Applies only if using '-w all' or '-w label_tick'
                Indicates a factor to decrease/increase the size of labels
                For example, to use a value 1.5 times bigger than normal use '-x 1.5', or to reduce them by half use '-x 0.5'
                Default = '1'"),

  make_option(c("-d", "--graphics_device"), default="pdf",
              help="Indicates if the plots should be either 'pdf' (default) or 'png'
                or type 'NA' to only generate text files with AUC values
                Default = 'pdf'")
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat            <- opt$infile_mat
InfileColours        <- opt$infile_colours
Outdir               <- opt$outdir
PrefixOutfiles       <- opt$prefix_outfiles
Lwd                  <- opt$lwd
PrintLabelsAndTicks  <- opt$print_plot_ticks_labels_legend
CexLabels            <- as.numeric(opt$cex_labels)
GraphicsDevice       <- opt$graphics_device

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
### Load data, get plots and AUC's
####################################
writeLines("\n*** Load data, get plots and AUC's ***\n")

mat<-read.table(InfileMat, header = T, row.names = 1, sep = "\t")

####################################
### Define colours
####################################
writeLines("\n*** Define colours ***\n")

DatasetLabelToColourHex  = list()

if (regexpr("^NA$", InfileColours, ignore.case = T)[1] == 1) {
  writeLines("\n*** Using default colours ***\n")
  
    ColoursNumbToHex <- list()
    ColoursNumbToHex <- list(
    ### The order of line colours in plots will follow the order of colours in list(ColoursNumbToHex)
    "#E69F00",    ## orange        
    "#009E73",    ## bluishgreen   
    "#CC79A7",    ## reddishpurple 
    "#1E90FF",    ## dodgerblue    
    "#D55E00",    ## vermillion    
    "#8B8989",    ## snow4         
    "#FFD700",    ## yellow        
    "#43CD80",    ## seagreen3     
    "#FFC0CB",    ## pink	      	
    "#56B4E9",    ## skyblue       
    "#8B4789",    ## orchid4       
    "#0072B2",    ## blue	      	
    "#000000"     ## black
  )
  
  for (d in 2:ncol(mat)) {
    DatasetLabelToColourHex[[colnames(mat)[d]]] <- as.character(ColoursNumbToHex[d-1])
  }
  
}else{
  writeLines("\n*** Using colours defined by --infile_colours ***\n")
  
  colours.mat<-read.table(InfileColours, header = F, row.names = 1, sep = "\t", comment.char = "")
  colnames(colours.mat)<-c("Colour_HEX")
  DatasetLabelToColourHex = list()
  for (d in 1:nrow(colours.mat)) {
    DatasetLabelToColourHex[[rownames(colours.mat)[d]]] <- as.character(colours.mat[d,"Colour_HEX"][[1]])
  }
}

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
  
  ### Get colour for plot 
  if ((colheaderPrediction %in% names(DatasetLabelToColourHex)) == TRUE) {
    AUCs_ROC_cols[i]<-DatasetLabelToColourHex[[colheaderPrediction]]
    AUCs_PR_cols[i] <-DatasetLabelToColourHex[[colheaderPrediction]]
  }else{
    stop(c("ERROR!!! couldn't get colour definition for column header", colheaderPrediction, "\n"))
  }
  
  ## Get ROC and PR 'saw-toothed' curve coordinates and AUC values
  curves.r <- evalmod(scores = predsToPlot, labels = labelsToPlot)
  aucs.r   <- auc(curves.r)
  AUCs_ROC_vals[i] <-round(aucs.r$aucs[[1]], digits = 2)
  AUCs_PR_vals[i]  <-round(aucs.r$aucs[[2]], digits = 2)
  AUCs_ROC_text[i]<-paste(colnames(mat)[[i]], "=", (round(aucs.r$aucs[[1]], digits = 2)), sep = "", collapse = "")
  AUCs_PR_text[i] <-paste(colnames(mat)[[i]], "=", (round(aucs.r$aucs[[2]], digits = 2)), sep = "", collapse = "")
  ### This forces lty to start from 2, because lty=1 is a solid line,
  ### which doesn't allow to see overlapping lines
  LineTypes[i] <- as.integer((i - 2) / length(names(DatasetLabelToColourHex))) + 2
  
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
      plot(x=curves.r$rocs[[1]]$x, y=curves.r$rocs[[1]]$y, lwd=Lwd, col=AUCs_ROC_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', xlab = "", ylab = "")
      
      if (grepl(pattern = "all|label_tick", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        mtext(text = "1 - Specificity", side=1, line = 3.8, cex=CexLabels)
        mtext(text = "Recall",          side=2, line = 2.2, cex=CexLabels)
        axis(side = 1, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels), mgp=c(3,2,0))
        axis(side = 2, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels))
        par(cex=1)
      }else if (grepl(pattern = "ticks", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        axis(side = 1, at=c(0,0.5,1), labels = FALSE)
        axis(side = 2, at=c(0,0.5,1), labels = FALSE)
      }

      ### Guide line
      par(new=T)
      plot(c(1:0), c(1:0), lty=1, col="gray60", lwd=1, type="l", axes=FALSE, bty="n", xlab="", ylab="", main="")
      
    }else{
  		par(new=T)
      plot(x=curves.r$rocs[[1]]$x, y=curves.r$rocs[[1]]$y, lwd=Lwd, col=AUCs_ROC_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }

    ## Plot PR 'saw-toothed' curves
    dev.set(3)
    if (i == 2) {
      plot(x=curves.r$prcs[[1]]$x, y=curves.r$prcs[[1]]$y, lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', xlab = "", ylab = "")
      
      if (grepl(pattern = "all|label_tick", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        mtext(text = "Recall",    side=1, line = 3.8, cex=CexLabels)
        mtext(text = "Precision", side=2, line = 2.2, cex=CexLabels)
        axis(side = 1, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels), mgp=c(3,2,0))
        axis(side = 2, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels))
        par(cex=1)
      }else if (grepl(pattern = "ticks", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        axis(side = 1, at=c(0,0.5,1), labels = FALSE)
        axis(side = 2, at=c(0,0.5,1), labels = FALSE)
      }
      
      ### Guide line
      abline(h=0.5, lty=1, col="gray60", lwd=1)
      
    }else{
      par(new=T)
      plot(x=curves.r$prcs[[1]]$x, y=curves.r$prcs[[1]]$y, lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }
    
    ## Plot PR 'interpolated' curves
    dev.set(4)
    if (i == 2) {
      plot(x=x.values[2:length(x.values)],y=cm[length(cm):1],lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', xlab = "", ylab = "")
      
      if (grepl(pattern = "all|label_tick", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        mtext(text = "Recall",    side=1, line = 3.8, cex=CexLabels)
        mtext(text = "Precision", side=2, line = 2.2, cex=CexLabels)
        axis(side = 1, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels), mgp=c(3,2,0))
        axis(side = 2, at=c(0,0.5,1), labels = c(0,"",1), par(cex=CexLabels))
        par(cex=1)
      }else if (grepl(pattern = "ticks", ignore.case = T, x = PrintLabelsAndTicks) == T) {
        axis(side = 1, at=c(0,0.5,1), labels = FALSE)
        axis(side = 2, at=c(0,0.5,1), labels = FALSE)
      }

      ### Guide line
      abline(h=0.5, lty=1, col="gray60", lwd=1)
      
    }else{
      par(new=T)
      plot(x=x.values[2:length(x.values)],y=cm[length(cm):1], lwd=Lwd, col=AUCs_PR_cols[i][[1]], lty=LineTypes[i][[1]], type="l", xlim=c(0,1), ylim=c(0,1), axes=FALSE,bty="n",xlab="",ylab="",main="")
    }
  }
}

if (GeneratePlots == 1) {
  
  if (grepl(pattern = "all|legend", ignore.case = T, x = PrintLabelsAndTicks) == T) {
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
  }
  
  ### Close graphic devices
  graphics.off()
}

####################################
### Generate AUC's outfiles
####################################
writeLines("\n*** Generate AUC's outfiles ***\n")

write(x=paste("Dataset", "ROC_AUC", sep = "\t", collapse = ""), file = OutRocTsv)
write(x=paste(colnames(mat)[2:length(colnames(mat))], AUCs_ROC_vals[2:length(colnames(mat))], sep = "\t", collapse = "\n"), append = T, file = OutRocTsv)

write(x=paste("Dataset", "PR_AUC", sep = "\t", collapse = ""), file = OutPrStTsv)
write(x=paste(colnames(mat)[2:length(colnames(mat))], AUCs_PR_vals[2:length(colnames(mat))], sep = "\t", collapse = "\n"), append = T, file = OutPrStTsv)

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

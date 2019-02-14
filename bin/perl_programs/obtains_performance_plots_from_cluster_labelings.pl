#!/usr/bin/perl

########################
### See documentation below, in section 'START INSTRUCTIONS TO RUN THIS PROGRAM'
########################

########################
### EXTERNAL DEPENDENCIES
########################
# R (Rscript)
# R script 'obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R' and its dependencies
# Perl
# Perl modules indicated below by 'use'
##########################

############################
#### Step 1 -- Define parameters, default modules, infiles, and other dependencies
############################

use LoadParameters::Parameters;
use ReformatPerlEntities::ObtainOutfileWOpath;
use PathsDefinition::PathsToInputs;
use ReformatNumbers::NAvalues;

### Change here the path to file 'obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R'
$PathToScriptForROCandPR = "~/r_programs/obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R";

$ThisProgramName = $0;
$ThisProgramName =~ s/\S+\///;

$CommentsForHelp = "
#####################################################################################
################### START INSTRUCTIONS TO RUN THIS PROGRAM ##########################
###
### Entering a -infile_list_infiles of cell cluster labeling results from methods like CIBERSORT, GSEA, GSVA, ORA, etc.
### and a -infile_gold_standards with reference cluster annotations, this script plots Receiver Operating Characteristic (ROC)
### and Precision / Sensitivity (PR) curves, and reports their Area Under the Curve (AUC)
###
### -------------------------------------------INPUTS----------------------------------------------
###
### [1]
### -infile_list_infiles is a list of outfiles from cell cluster labeling methods, dataset IDs, and infile types, like:
###
### #Outfile from cell cluster labeling method       Dataset ID    Infile type
### #----------------------------------------------  ------------  -----------
### Infile_1*CIBERSORT_enrichment_scores.tsv         Unique_ID_1   Generic_ES
### Infile_2*Pval.Unfiltered.mat.txt                 Unique_ID_2   GSEA_Pvalues
### Infile_3*GSVA_all_scores_table.tsv               Unique_ID_3   GSVA_ES
### Infile_4*ORA.PvalUncorrected.All.mat             Unique_ID_4   ORA_Pvalues
###
### Infile types accepted                        Presumably produced by script
### ------------------------------------------   ---------------------------------------
### a) *CIBERSORT_enrichment_scores.tsv          obtains_CIBERSORT_for_MatrixColumns.pl
### b) *Pval.Unfiltered.mat.txt                  obtains_GSEA_for_MatrixColumns.pl
### c) *GSVA_all_scores_table.tsv                obtains_GSVA_for_MatrixColumns.R
### d) *ORA.PvalUncorrected.All.mat              obtains_ORA_for_MatrixColumns.pl
###
### Note:
### Infile type 'Generic_ES' option can be used to process any matrices of clusters (rows) vs. cell types (columns), like:
###
### Cluster  B_CELLS  T_CELLS  HEPATOCYTE
### clust1   0.05     0.01     0.9
### clust2   0.05     0.02     0.8
### clust3   0.01     0.9      0.01
###
### And hence could be used to process other cell cluster labeling methods not evaluated by Diaz-Mejia JJ et al (2019)
###
###
### [2]
### -infile_gold_standards in format like:
###
### #Cluster  Cell type
### #-------  ---------
### clust1    HEPATOCYTE
### clust2    HEPATOCYTE
### clust3    T_CELLS
###
### ----------------------------------------MAIN OUTPUTS-------------------------------------------
###
### [1]
### *merged.tsv table with compiled scores from inputted method results and gold standards
###
### [2]
### PERFORMANCE_PLOTS/*.auc table with Area Under the Curve values
###
### [2]
### PERFORMANCE_PLOTS/*.pdf table with the ROC and PR curve plots
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles         (path/name to the directory where outfiles will be saved)
###   -infile_list_infiles   (path/name to the list of infiles to process)
###   -infile_gold_standards (path/name to the gold standards)
###   -use_graphics_device   (indicates if the type of outfile: 'pdf' or 'png
###                           Or type 'NA' to get only the Area Under the Curve (AUC) values, without plots)
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

&Readme;
&Parameters;

$useForLogOf0 = 0.0000000000001; ### To be use for -log(0) of p-values

############################
#### Step 2 -- Loading data
############################

### Load cluster labeling results

$infile_list = $hashParameters{infile_list_infiles};
$infile_list =~ s/~\//\/$Users_home\/$DefaultUserName\//;

if ($infile_list =~ /\.bz2$/) {
open (INFILE_LIST, "bzip2 -qdc $infile_list |") or die "Can't open '$infile_list' (expected bzipped infile_list)\n";
}elsif ($infile_list =~ /\.gz$/) {
open (INFILE_LIST, "gzip -qdc $infile_list |") or die "Can't open '$infile_list' (expected gzipped infile_list)\n";
}else{
open (INFILE_LIST, "<$infile_list") or die "Can't open '$infile_list' (expected unzipped infile_list)\n";
}

$countDatasets = 0;
while ($line = <INFILE_LIST>) {
chomp $line;

	unless ($line =~ /^#/) {
		if ($line =~ /^(\S+)(\s+)(\S+)(\s+)(\S+)$/) {
		$infile_clust_labels   = $1;
		$DatasetID             = $3;
		$infile_type           = $5;
		$infile_clust_labels   =~ s/~\//\/$Users_home\/$DefaultUserName\//;
		$infile_gold_standards =~ s/~\//\/$Users_home\/$DefaultUserName\//;
		}else{
		die "\n\nERROR!!! unexpected format in line:\n'$line'\n\n";
		}
		
		if ($hashAllDatasetIDs{$DatasetID}) {
		die "\n\nERROR!!! 'unique_flag_dataset' strings must be unique, but '$DatasetID' appears more than once in -infile_list_infiles\n\n";
		}else{
		$countDatasets++;
		$hashAllOrderedDatasets{$countDatasets} = $DatasetID;
		$hashAllDatasetIDs{$DatasetID} = 1;
		}
	
		if (-f $infile_clust_labels) {
			if ($infile_type =~ /^GSVA_ES$/i) {
			&LoadDataForBenchmark_GSVA($infile_clust_labels,$DatasetID);
			
			}elsif ($infile_type =~ /^ORA_Pvalues$/i) {
			&LoadDataForBenchmark_ORA($infile_clust_labels,$DatasetID);

			}elsif ($infile_type =~ /^GSEA_Pvalues$/i) {
			&LoadDataForBenchmark_GSEA($infile_clust_labels,$DatasetID);
			
			}elsif ($infile_type =~ /^Generic_ES$/i) {
			&LoadDataForBenchmark_GENERIC_ES($infile_clust_labels,$DatasetID);
			
			}else{
			die "\n\nERROR!!! unexpected infile_type '$infile_type'\n\n";
			}
		}
	}
}
close INFILE_LIST;


### Load gold standards

$infile_gold_standards = $hashParameters{infile_gold_standards};
$infile_gold_standards =~ s/~\//\/$Users_home\/$DefaultUserName\//;
&LoadDataForBenchmark_Golds($infile_gold_standards,$DatasetID);


############################
#### Step 3 -- Generate outfile
############################

open OUTFILE, ">$hashParameters{path_outfiles}/$outfileWOpath.merged.tsv" or die "Can't open '$hashParameters{path_outfiles}/$outfileWOpath.merged.tsv'\n";

print OUTFILE "ClusterId_ClassId\tLabel";
$hashAllOrderedDatasets{$countDatasets} = $DatasetID;

foreach $NumberOfDataset (1..$countDatasets) {
$DatasetID = $hashAllOrderedDatasets{$NumberOfDataset};
	foreach $ColumnToCompare (sort keys %{$hashAllScoresPerDataset{$DatasetID}}) {
	print OUTFILE "\t$DatasetID" . "__" . "$ColumnToCompare";
	}
}
print OUTFILE "\n";

foreach $clusterid (sort keys %hashAllClusterIds) {
	foreach $classid (sort keys %hashAllClasses) {
		
		if ($hashGoldStandards{$clusterid}{$classid}) {
		$Label = 1;
		}else{
		$Label = 0;
		}
		print OUTFILE "$clusterid.vs.$classid\t$Label";

		foreach $NumberOfDataset (1..$countDatasets) {
		$DatasetID = $hashAllOrderedDatasets{$NumberOfDataset};
			foreach $ColumnToCompare (sort keys %{$hashAllScoresPerDataset{$DatasetID}}) {
				if ($hashAllScoresPerDataset{$DatasetID}{$ColumnToCompare}{$clusterid}{$classid}) {
				$score = $hashAllScoresPerDataset{$DatasetID}{$ColumnToCompare}{$clusterid}{$classid};
				$score =~ s/^_//;
				}else{
				$score = "NA";
				}
			print OUTFILE "\t$score";
			}
		}
		print OUTFILE "\n";
	}
}
close OUTFILE;

############################
#### Step 4 -- Generate plots
############################

system "Rscript $PathToScriptForROCandPR -i $hashParameters{path_outfiles}/$outfileWOpath.merged.tsv -o $hashParameters{path_outfiles}/PERFORMANCE_PLOTS -p $outfileWOpath -d $hashParameters{use_graphics_device} -l 5";

&PrintParameters;

print "\n\n  Done!!!\n  Check '$hashParameters{path_outfiles}/$outfileWOpath.*' for outfiles\n\n";

exit;

########################################################
################ END OF PROGRAM ########################
########################################################


########################################################
################ START SUBROUTINES #####################
########################################################


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Parameters {

########## Print "Usage" for user

print "$CommentsForHelp\n\n";

##########################
######## Options and Infiles

use Cwd 'abs_path';
$ScriptName = abs_path($0);
$Parameters .= "$ScriptName\n";

chomp @ARGV;
@arrayInputtedOneLineCommands = @ARGV;

%hashParametersTolookFor = (
'path_outfiles' => 1,
'infile_list_infiles' => 1,
'infile_gold_standards' => 1,
'use_graphics_device' => 1,
);

#######################
#### Starts -- Evaluate parameters

&LoadParameters::Parameters::MainSubParameters(\%hashParametersTolookFor,\@arrayInputtedOneLineCommands);
$Parameters .= "$MoreParameters";

## Defining prefix string for OUTFILE
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_list_infiles});

#### Ends -- Evaluate parameters
#######################

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub PrintParameters {

### Printing out parameters. Need to be concatenated at sub(Parameters)
open PARAMETERS, ">$hashParameters{path_outfiles}/$outfileWOpath.PlotClusterLabelingPerformance.Parameters" or die "Can't open '$hashParameters{path_outfiles}/$outfileWOpath.PlotClusterLabelingPerformance.Parameters'\n";
print PARAMETERS "$Parameters";
close PARAMETERS;
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Readme {

my ($date) = `date`;
chomp $date;

$Parameters .= "
#################################################################
# Javier Diaz -- $date
# javier.diazmejia\@gmail.com
#################################################################\n
$Extras
#################################################################
######################### PARAMETERS ############################
#################################################################\n\n";

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadDataForBenchmark_GSVA {

($infile,$DatasetID) = @_;

print "Loading '$DatasetID' dataset from:\n$infile\n";

open (INFILE, "<$infile") or die "Can't open '$infile' infile by LoadDataForBenchmark_GSVA\n";

$countRows = 0;
while ($line = <INFILE>) {
chomp $line;
$c = 0;
	unless ($line =~ /^#/) {
	$countRows++;
		if ($countRows == 1) {
		@ColumnsHeaders = split ("\t", $line);
		$cMinusOne = -1;
			foreach $columnheader (@ColumnsHeaders) {
			$cMinusOne++;
				if ($cMinusOne == 0) {
				$cMinusOne = "0.0";
				}
			$hashColumnNumbersMinusOne{$DatasetID}{$columnheader} = $cMinusOne;
			}
		}else{
		@Data = split ("\t", $line);
		$classid   = @Data[$hashColumnNumbersMinusOne{$DatasetID}{CLASS}];
		$clusterid = @Data[$hashColumnNumbersMinusOne{$DatasetID}{ColumnHeader}];
		$es        = @Data[$hashColumnNumbersMinusOne{$DatasetID}{EnrichmentScore}];
		$hashAllScoresPerDataset{$DatasetID}{GSVA_ES}{$clusterid}{$classid} = "_$es";
		$hashAllClusterIds{$clusterid} = 1;
		$hashAllClasses{$classid} = 1;
		}
	}
}
close INFILE;

}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadDataForBenchmark_ORA {

($infile,$DatasetID) = @_;

print "Loading '$DatasetID' dataset from:\n$infile\n";

open (INFILE, "<$infile") or die "Can't open '$infile' infile by LoadDataForBenchmark_ORA\n";

$countRows = 0;
while ($line = <INFILE>) {
chomp $line;
$c = 0;
	unless ($line =~ /^#/) {
	$countRows++;
		if ($countRows == 1) {
		@ColumnsHeaders = split ("\t", $line);
		$cMinusOne = -1;
			foreach $columnheader (@ColumnsHeaders) {
			$cMinusOne++;
				if ($cMinusOne == 0) {
				$cMinusOne = "0.0";
				}
			$hashColumnNumbersMinusOne{$DatasetID}{$columnheader} = $cMinusOne;
			$hashColumnHeadersMinusOne{$DatasetID}{$cMinusOne} = $columnheader;
			}
		}else{
		@Data = split ("\t", $line);
		$cMinusOne = -1;
		$enter_this_row = 0;

			foreach $d (@Data) {
			$cMinusOne++;
			$columnheader = $hashColumnHeadersMinusOne{$DatasetID}{$cMinusOne};
				
				if ($cMinusOne == 0) {
				$classid = @Data[$cMinusOne];
					
					if ($classid =~ /^(cutoff_pos)(---)(\S+)(---)(\S+)$/) {
					$classid = $3;
					$enter_this_row = 1;
					}
						
				}else{
					if ($enter_this_row == 1) {
						unless ($hashNAvalues{$d}) {
						$clusterid = $columnheader;
							if ($d == 0) {
							$d = $useForLogOf0;
							}
						$mLogPvalue = log($d) * -1;
						$hashAllScoresPerDataset{$DatasetID}{ORA_mLogPval}{$clusterid}{$classid} = "_$mLogPvalue";
						$hashAllClusterIds{$clusterid} = 1;
						$hashAllClasses{$classid} = 1;
						}
					}
				}
			}
		}
	}
}
close INFILE;

}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadDataForBenchmark_GSEA {

($infile,$DatasetID) = @_;

print "Loading '$DatasetID' dataset from:\n$infile\n";

open (INFILE, "<$infile") or die "Can't open '$infile' infile by LoadDataForBenchmark_GSEA\n";

$countRows = 0;
while ($line = <INFILE>) {
chomp $line;
$c = 0;
	unless ($line =~ /^#/) {
	$countRows++;
		if ($countRows == 1) {
		@ColumnsHeaders = split ("\t", $line);
		$cMinusOne = -1;
			foreach $columnheader (@ColumnsHeaders) {
			$cMinusOne++;
				if ($cMinusOne == 0) {
				$cMinusOne = "0.0";
				}
			$hashColumnNumbersMinusOne{$DatasetID}{$columnheader} = $cMinusOne;
			$hashColumnHeadersMinusOne{$DatasetID}{$cMinusOne} = $columnheader;
			}
		}else{
		@Data = split ("\t", $line);
		$cMinusOne = -1;
		$enter_this_row = 0;

			foreach $d (@Data) {
			$cMinusOne++;
			$columnheader = $hashColumnHeadersMinusOne{$DatasetID}{$cMinusOne};
				
				if ($cMinusOne == 0) {
				$classid = @Data[$cMinusOne];
					
					if ($classid =~ /^(pos)(\.)(\S+)(\.\d+\.)(\S+)$/) {
					$classid = $5;
					$enter_this_row = 1;
					}
						
				}else{
					if ($enter_this_row == 1) {
						unless ($hashNAvalues{$d}) {
						$clusterid = $columnheader;
						
							if ($d == 0) {
							$d = $useForLogOf0;
							}
						
						$mLogPvalue = log($d) * -1;
						$hashAllScoresPerDataset{$DatasetID}{GSEA_mLogPval}{$clusterid}{$classid} = "_$mLogPvalue";
						$hashAllClusterIds{$clusterid} = 1;
						$hashAllClasses{$classid} = 1;
						}
					}
				}
			}
		}
	}
}
close INFILE;

}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadDataForBenchmark_GENERIC_ES {

($infile,$DatasetID) = @_;

print "Loading '$DatasetID' dataset from:\n$infile\n";

open (INFILE, "<$infile") or die "Can't open '$infile' infile by LoadDataForBenchmark_GENERIC_ES\n";

$countRows = 0;
while ($line = <INFILE>) {
chomp $line;
$c = 0;
	unless ($line =~ /^#/) {
	$countRows++;
		if ($countRows == 1) {
		@ColumnsHeaders = split ("\t", $line);
		$cMinusOne = -1;
			foreach $columnheader (@ColumnsHeaders) {
			$columnheader =~ s/ /_/g;
			$cMinusOne++;
				unless ($cMinusOne == 0) {
				$hashColumnNumbersMinusOne{$DatasetID}{$cMinusOne} = $columnheader;
				}
			}
		}else{
		@Data = split ("\t", $line);
		$cMinusOne = -1;
		$clusterid = @Data[0];
		
			foreach $d (@Data) {
			$cMinusOne++;
				unless ($cMinusOne == 0) {
				$classid   = $hashColumnNumbersMinusOne{$DatasetID}{$cMinusOne};
				$es        = $d;
				$hashAllScoresPerDataset{$DatasetID}{Generic_ES}{$clusterid}{$classid} = "_$es";
				$hashAllClusterIds{$clusterid} = 1;
				$hashAllClasses{$classid} = 1;
				}
			}
		}
	}
}
close INFILE;

}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadDataForBenchmark_Golds {

($infile) = @_;

print "Loading gold standards from:\n$infile\n";

open (INFILE, "<$infile") or die "Can't open '$infile' infile by LoadDataForBenchmark_Golds\n";

while ($line = <INFILE>) {
chomp $line;
	unless ($line =~ /^#/) {
		if ($line =~ /^(\S+)(\s+)(\S+)/) {
		$clusterid = $1;
		$classid   = $3;
		
		## Indexed in this way in principle should allow to have more than one
		## gold standard 'class' annotation per cluster but I'm using only one anyways
		$hashGoldStandards{$clusterid}{$classid} = 1;
		}
	}
}
close INFILE;

}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#!/usr/bin/perl

########################
### See documentation below, in section 'START INSTRUCTIONS TO RUN THIS PROGRAM'
########################

########################
### EXTERNAL DEPENDENCIES
########################
### (1)
### gsea-3.0.jar  -- (GSEA standalone executable java *.jar file)
### It can be downloaded from http://software.broadinstitute.org/gsea/downloads.jsp
### For one-line commands type:
### java -cp ~/.../gsea-3.0.jar -Xmx512m xtools.gsea.GseaPreranked  -help
###
### In this script specify the path/name to the gsea-3.0.jar by variable $GseaExecutable
###
### (3)
### Java (tested with version 1.8.0_162)
###
### (4)
### Perl modules -- Defined by 'use' command below
########################

############################
#### Step 1 -- Define parameters, default modules, infiles, and other dependencies
############################

use LoadParameters::Parameters;
use ReformatPerlEntities::ObtainOutfileWOpath;
use PathsDefinition::PathsToInputs;
use LoadData::LoadClasses;
use LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue;
use LoadData::LoadMatrixNoStrain_FlagDataset_Undirected;
use Date::Calc qw(Delta_DHMS);

### Change here the path to file 'gsea-3.0.jar'
$GseaExecutable = "~/PROGRAMS/GSEA/gsea-3.0.jar";

$ThisProgramName = $0;
$ThisProgramName =~ s/\S+\///;

$CommentsForHelp = "
#####################################################################################
################### START INSTRUCTIONS TO RUN THIS PROGRAM ##########################
###
### Runs GSEA for each column in -infile_matrix of genes (rows) vs. arrays (columns) e.g. samples, clusters, etc
### using an -infile_classes with gene-classes (e.g. cell-type markers, Gene Ontology terms, etc)
###
### To compute enrichment, the user can opt to use either the INTERSECTION of -infile_matrix and -infile_classes
### or ALL genes in -infile_matrix as the universe of genes
###
### -------------------------------------------INPUTS----------------------------------------------
###
### [1]
### a -infile_matrix in <tab> delimited format, like:
### COMPLETE   Array_1   Array_2   Array_3
### Gene1      0.01      0.01      0.03
### Gene2     -0.01      NODATA   -0.15
### Gene3     NODATA     0.05      0.7
### ...etc
###
### [2]
### a -infile_classes in <tab> delimited *.gmt format, like:
### CellType_1_id  CellType_1_name  Gene1  Gene2  Gene3
### CellType_2_id  CellType_2_name  Gene4  Gene5
### CellType_3_id  CellType_3_name  Gene5  Gene6
### ...etc
###
### ----------------------------------------MAIN OUTPUTS-------------------------------------------
###
### [1]
### Tables of P-value and FDR values obtained from GSEA for each column of -infile_matrix vs. each gene-class of -infile_classes
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles      (path/name to the directory where outfiles will be saved)
###   -prefix_outfiles    (a string for the outfile names)
###   -infile_matrix      (path/name to the matrix of genes vs. conditions)
###   -infile_classes     (path/name to either the gene classes)
###   -nperm              (number of permutations to generate GSEA statistical analyses. GSEA developers recommend minimum '1000')
###   -use_universe       (indicates set of genes to use as universe from -infile_matrix and/or -infile_classes
###                        Either 'n1' to use all genes in -infile_matrix,
###                        Or 'i' to use the intersection of -infile_matrix and -infile_classes)
###   -cutoff_print_p     (once GSEA is computed, this cutoff selects Nom P-values that will be printed out in summary outfiles (e.g. 0.01))
###                        Note: also will generate outfiles without filtering neither P-values nor Q-values
###   -cutoff_print_q     (once GSEA is computed, this cutoff selects FDR Q-values that will be printed out in summary outfiles (e.g. 0.25))
###                        Note: also will generate outfiles without filtering neither P-values nor Q-values
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

### This will be used to get overall computing time
($year_overall1,$month_overall1,$day_overall1,$hour_overall1,$minute_overall1,$second_overall1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

&Readme;
&Parameters;

$GseaExecutable =~ s/~\//\/$Users_home\/$DefaultUserName\//;

unless (-f $GseaExecutable) {
die "\n\nERROR!!! Couldn't find gsea*jar executable at '$GseaExecutable'\n"
}

## By default GSEA plots only the top 20 sets (including detailed reports on 'core' gene sets, necessary to generate gene-level reports
## This may not be enough to gather all info needed for the full report on all enriched sets
## thus we'll first plot a few sets ($PlotTopFew) and if they are not enough then plot many ($PlotTopMany)
## In any case this shouldn't affect the class-level enrichment outfiles
$PlotTopFew  = 50;  
$PlotTopMany = 2000;

unless ($hashParameters{use_universe} =~ /^(i|n1)$/i) {
die "\n\nERROR!!! unexpected option '-use_universe $hashParameters{use_universe}'\n\n";
}

############################
#### Step 1 -- Loading data and preparing infiles for GSEA
############################

### NOTE: load -infile_classes before -infile_matrix in case '-use_universe i' is used

LoadData::LoadClasses::LoadClasses($hashParameters{infile_classes},"ALL",1,1000000,"NA","NA","yes");
LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue::LoadMatrixReturnRnkEachColumnDepurateByMaxValue($hashParameters{path_outfiles},$hashParameters{infile_matrix},"NA","NA","NA",$hashParameters{use_universe},\%HashEachKeyClassesPassingCutoff);

############################
#### Step 2 -- Sending each column to GSEA
############################

#Removing previous (presumably aborted) results
$outdir = "$hashParameters{path_outfiles}/GSEA/$hashParameters{prefix_outfiles}";

if (-d $outdir) {
`rm -r $outdir`;
print "\n\nWARNING!!! removing preexisting outdir:\n'$outdir'\n";
}
system "mkdir -p $outdir";

print "Performing GSEA for each column\n";

### This will be used to get GSEA computing time
($year_gsea1,$month_gsea1,$day_gsea1,$hour_gsea1,$minute_gsea1,$second_gsea1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

foreach $c (1..$NumberOfColumnHeaders) {
$ColumnName = $hashColumnsNames{$c};
print "\nPerforming '$ColumnName' of '$hashParameters{infile_matrix}' vs. $hashParameters{infile_classes} ($ColumnName.rnk)'\n";
$InFileName  = "$hashParameters{path_outfiles}/$ColumnName.rnk";
$OutFileName = "$ColumnName.vs.$outfileWOpath_classes.GSEA";
&SendToGseaInBatchVerifyingOutfilesForGeneLevel($InFileName,$OutFileName,$PlotTopFew);
}

### This will be used to get GSEA computing time
($year_gsea2,$month_gsea2,$day_gsea2,$hour_gsea2,$minute_gsea2,$second_gsea2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

############################
#### Step 3 -- Printing summary out
############################

print "Printing summary\n";

open SUMMARYPVALALL, ">$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Pval.Unfiltered.mat.txt"                                                        or die "Can't open \"$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Pval.Unfiltered.mat.txt\"\n";
open SUMMARYPVALFIL, ">$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Pval.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.mat.txt" or die "Can't open \"$outdir/ClassLevel/$hashParameters{prefix_outfiles}.GSEA.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.Pval.mat.txt\"\n";
open SUMMARYFDRALL,  ">$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Qval.Unfiltered.mat.txt"                                                        or die "Can't open \"$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Qval.Unfiltered.mat.txt\"\n";
open SUMMARYFDRFIL,  ">$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Qval.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.mat.txt" or die "Can't open \"$outdir/ClassLevel/$hashParameters{prefix_outfiles}.Qval.GSEA.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.mat.txt\"\n";

print SUMMARYPVALALL "GSEA.NomPvalue.All";
print SUMMARYPVALFIL "GSEA.NomPvalue.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}";
print SUMMARYFDRALL  "GSEA.FDRQvalue.All";
print SUMMARYFDRFIL  "GSEA.FDRQvalue.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}";

foreach $c (1..$NumberOfColumnHeaders) {
	if ($hashColumnsNames{$c}) {
	$ColumnName = $hashColumnsNames{$c};
	$NumberOfGenes = 0;
		if ($hashCountScoresPassingCutoffsWONas{$ColumnName}) {
		$NumberOfGenes = $hashCountScoresPassingCutoffsWONas{$ColumnName};
		}else{
		$NumberOfGenes = "NA";
		}
	
	print SUMMARYPVALALL "\t$ColumnName";
	print SUMMARYPVALFIL "\t$ColumnName";
	print SUMMARYFDRALL  "\t$ColumnName";
	print SUMMARYFDRFIL  "\t$ColumnName";
	}else{
	die "\n\nERROR!!! couldn't find $hashColumnsNames{$c}\n";
	}
}
print SUMMARYPVALALL "\n";
print SUMMARYPVALFIL "\n";
print SUMMARYFDRALL  "\n";
print SUMMARYFDRFIL  "\n";
	
foreach $TypeOfScoreClass (sort keys %hashAllClassesTypes) {
($typeOfScore,$classid) = split ("\t", $TypeOfScoreClass);
$ClassDetailsPartA = "";
$ClassDetailsPartB = "";
$ToPrintConcatenatedCpvalAll = "";
$ToPrintConcatenatedCqvalAll = "";
$ToPrintConcatenatedCpvalFiltered = "";
$ToPrintConcatenatedCqvalFiltered = "";

	if ($HashClassRename{$classid}) {
	$ClassName = $HashClassRename{$classid};
	}else{
	$ClassName = $ClassName;
	}

$ClassDetailsPartA = "$typeOfScore.$classid";

	foreach $c (1..$NumberOfColumnHeaders) {
		if ($hashColumnsNames{$c}) {
		$ColumnName = $hashColumnsNames{$c};

		$PassPvalueFilter = 0;
		$PassQvalueFilter = 0;
		
			if ($hashSummary{$ColumnName}{$TypeOfScoreClass}{nompval}) {
			$nompval = $hashSummary{$ColumnName}{$TypeOfScoreClass}{nompval};
			$nompval =~ s/_//;
				if ($hashSummary{$ColumnName}{$TypeOfScoreClass}{real}) {
				$size = $hashSummary{$ColumnName}{$TypeOfScoreClass}{real};
				}
				if ($nompval <= $hashParameters{cutoff_print_p}) {
				$PassPvalueFilter = 1;
				}
			}else{
			$nompval = "NA";
			}
	
			if ($hashSummary{$ColumnName}{$TypeOfScoreClass}{fdrqval}) {
			$fdrqval = $hashSummary{$ColumnName}{$TypeOfScoreClass}{fdrqval};
			$fdrqval =~ s/_//;
				if ($fdrqval <= $hashParameters{cutoff_print_q}) {
				$PassQvalueFilter = 1;
				}
			}else{
			$fdrqval = "NA";
			}
	
			if (($PassPvalueFilter == 1) && ($PassQvalueFilter == 1)) {
			$ToPrintConcatenatedCpvalFiltered .= "\t$nompval";
			$ToPrintConcatenatedCqvalFiltered .= "\t$fdrqval";
			}elsif ($nompval =~ /\d/) {
			$ToPrintConcatenatedCpvalFiltered .= "\tNS";
			$ToPrintConcatenatedCqvalFiltered .= "\tNS";
			}else{
			$ToPrintConcatenatedCpvalFiltered .= "\t$nompval";
			$ToPrintConcatenatedCqvalFiltered .= "\t$fdrqval";
			}
		$ToPrintConcatenatedCpvalAll .= "\t$nompval";
		$ToPrintConcatenatedCqvalAll .= "\t$fdrqval";
		}
	}
$ClassDetailsPartB = "$size.$ClassName";
	unless ($hashYaPrintClassName{$classid}) {
	$hashYaPrintClassName{$classid} = 1;
	}
	
$ToPrintConcatenatedCpvalAll =~ s/^\t//;
$ToPrintConcatenatedCqvalAll =~ s/^\t//;
$ToPrintConcatenatedCpvalFiltered =~ s/^\t//;
$ToPrintConcatenatedCqvalFiltered =~ s/^\t//;

print SUMMARYPVALALL "$ClassDetailsPartA.$ClassDetailsPartB\t$ToPrintConcatenatedCpvalAll\n";
print SUMMARYFDRALL  "$ClassDetailsPartA.$ClassDetailsPartB\t$ToPrintConcatenatedCqvalAll\n";
print SUMMARYPVALFIL "$ClassDetailsPartA.$ClassDetailsPartB\t$ToPrintConcatenatedCpvalFiltered\n";
print SUMMARYFDRFIL  "$ClassDetailsPartA.$ClassDetailsPartB\t$ToPrintConcatenatedCqvalFiltered\n";
}
close SUMMARYPVALALL;
close SUMMARYFDRALL;
close SUMMARYPVALFIL;
close SUMMARYFDRFIL;

$lookOrignalColumnNamesFile = "$outfileWOpath_matrix.OriginalColNames.txt";
if (-f $lookOrignalColumnNamesFile) {
`mv $lookOrignalColumnNamesFile $hashParameters{path_outfiles}`;
}

############################
#### Step 4 -- Generating log file with parameters used in the run
############################

&PrintParameters;

############################
#### Step 5 -- Remove temporary files
############################

print "Removing temporary files\n";

`rm $hashParameters{infile_classes}.Filtered.gmt`;

$date = `date`;
if ($date =~ /^(\S+\s+)(\S+)(\s+)(\S+)/) {
$d1 = $2;
$d2 = $4;
$deletedate = "";
$add = 0;

	if ($d1 =~ /[a-z][a-z][a-z]/i) {
		if ($d2 =~ /^\d$/) {
		$deletedate = "$d1$add$d2";
		}else{
		$deletedate = "$d1$d2";
		}
	}elsif ($d2 =~ /[a-z][a-z][a-z]/i) {
		if ($d1 =~ /^\d$/) {
		$deletedate = "$d2$add$d1";
		}else{
		$deletedate = "$d2$d1";
		}
	}else{
	print "WARNING!!! Couldn't get date in format like 'May10' to remove temporary directory. Date '$date' was processed\n\n";
	}
$deletedate =~ tr/[A-Z]/[a-z]/;
`rmdir $deletedate`;
}

foreach $c (1..$NumberOfColumnHeaders) {
$ColumnName = $hashColumnsNames{$c};
$InFileName  = "$hashParameters{path_outfiles}/$ColumnName.rnk";
$OutFileName = "$ColumnName.vs.$outfileWOpath_classes.GSEA";
system "rm $InFileName";
system "rm -r $outdir/ClassLevel/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/";
}


############################
#### Step 7 -- Report run time
############################

### This will be used to get overall computing time
($year_overall2,$month_overall2,$day_overall2,$hour_overall2,$minute_overall2,$second_overall2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

### overall runtime
($days_overall_diff, $hours_overall_diff, $minutes_overall_diff, $seconds_overall_diff) =
Delta_DHMS($year_overall1,$month_overall1,$day_overall1,$hour_overall1,$minute_overall1,$second_overall1,
	   $year_overall2,$month_overall2,$day_overall2,$hour_overall2,$minute_overall2,$second_overall2);
	   
$day_overall_seconds     = $days_overall_diff    * 86400;
$hour_overall_seconds    = $hours_overall_diff   * 3600;
$minute_overall_seconds  = $minutes_overall_diff * 60;
$all_overall_seconds     = $day_overall_seconds + $hour_overall_seconds + $minute_overall_seconds + $seconds_overall_diff;

### gsea runtime
($days_gsea_diff, $hours_gsea_diff, $minutes_gsea_diff, $seconds_gsea_diff) =
Delta_DHMS($year_gsea1,$month_gsea1,$day_gsea1,$hour_gsea1,$minute_gsea1,$second_gsea1,
	   $year_gsea2,$month_gsea2,$day_gsea2,$hour_gsea2,$minute_gsea2,$second_gsea2);
	   
$day_gsea_seconds     = $days_gsea_diff    * 86400;
$hour_gsea_seconds    = $hours_gsea_diff   * 3600;
$minute_gsea_seconds  = $minutes_gsea_diff * 60;
$all_gsea_seconds     = $day_gsea_seconds + $hour_gsea_seconds + $minute_gsea_seconds + $seconds_gsea_diff;

open OUTRUNTIME, ">$outdir/$hashParameters{prefix_outfiles}.GSEA.Runtime" or die "Can't open '$outdir/$hashParameters{prefix_outfiles}.GSEA.Runtime'\n";
print OUTRUNTIME "gsea\t$all_gsea_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";
close OUTRUNTIME;

print "Took time:
gsea\t$all_gsea_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";

############################
#### Step 8 -- Finishing program
############################

print "\n  Done!!!\n\nGSEA conducted for:\n'$hashParameters{infile_matrix}'\nvs.\n'$hashParameters{infile_classes}'\n\nCheck directory:\n$outdir\n\n";

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
'infile_matrix' => 1,
'prefix_outfiles' => 1,
'cutoff_print_p' => 1,
'cutoff_print_q' => 1,
'infile_classes' => 1,
'nperm' => 1,
'use_universe' => 1,
);

#######################
#### Starts -- Evaluate parameters

&LoadParameters::Parameters::MainSubParameters(\%hashParametersTolookFor,\@arrayInputtedOneLineCommands);
$Parameters .= "$MoreParameters";

## Create directory for outifles
$hashParameters{path_outfiles} =~ s/\/$//;
unless (-d $hashParameters{path_outfiles}) {
`mkdir -p $hashParameters{path_outfiles}`;
}

## Defining prefix string for OUTFILE
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_matrix});
$outfileWOpath_matrix = $outfileWOpath;

if ($hashParameters{infile_classes} =~ /-/) {
die "\n\nERROR!!! GSEA can't handle en dashes '-' in file names:\n'$hashParameters{infile_classes}'\n\n";
}
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_classes});
$outfileWOpath_classes = $outfileWOpath;

#### Ends -- Evaluate parameters
#######################

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub PrintParameters {
### Printing out parameters. Need to be concatenated at sub(Parameters)
open PARAMETERS, ">$outdir/$hashParameters{prefix_outfiles}.GSEA.Parameters" or die "Can't open '$outdir/$hashParameters{prefix_outfiles}.GSEA.Parameters'\n";
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
sub SendToGseaInBatchVerifyingOutfilesForGeneLevel {

##################################
#Sending one-line-commands to GSEA
##################################

my ($InFileName,$OutFileName,$PlotTop) = @_;

### Here removing files from previous runs (if any)
$DirectoryFromPreviousRunName = "";
$DirectoryFromPreviousRunLook = "$outdir/ClassLevel/$OutFileName\.GseaPreranked\.[0-9]*[0-9]";

if (-d $DirectoryFromPreviousRunLook) {
$DirectoriesFromPreviousRunName = `ls -d $DirectoryFromPreviousRunLook`;
}

if ($DirectoriesFromPreviousRunName =~ /\S/) {
@DirectoriesFromPreviousRunName = split ("\n", $DirectoriesFromPreviousRunName);
	foreach $dir (@DirectoriesFromPreviousRunName) {
	system "rm -r $dir";
	}
}

$CommandToSend = "java -cp $GseaExecutable -Xmx1024m xtools.gsea.GseaPreranked \\
-gmx $hashParameters{infile_classes}.Filtered.gmt \\
-norm meandiv \\
-nperm $hashParameters{nperm} \\
-rnk $InFileName \\
-scoring_scheme weighted \\
-rpt_label $OutFileName \\
-make_sets true \\
-plot_top_x $PlotTop \\
-rnd_seed timestamp \\
-set_max 1000000 \\
-set_min 1 \\
-zip_report false \\
-out $outdir/ClassLevel \\
-gui false";

### These commands from gsea2.jar dissapear in gsea-3.0.jar
if ($GseaExecutable =~ /gsea2\.jar/) {
$CommandToSend .= " \\
-collapse false \\
-mode Max_probe \\
-include_only_symbols true";
}elsif ($GseaExecutable =~ /gsea-3\.0\.jar/) {
## Do nothing
}else{
die "\n\nERROR!!! unexpected gsea*jar executable from: '$GseaExecutable'\n\n";
}

print "\n\nRUNNING\n$CommandToSend\n\n";

system "$CommandToSend";

### Compiling results for summary
$baseForNegClassLevelReports = "";
$baseForPosClassLevelReports = "";
	
$resultsNeg = "";
$resultsPos = "";
$resultsNeg = `ls $outdir/ClassLevel/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/gsea_report_for_na_neg_[0-9]*[0-9]\.xls`;
$resultsPos = `ls $outdir/ClassLevel/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/gsea_report_for_na_pos_[0-9]*[0-9]\.xls`;

	if ($resultsNeg =~ /\S/ && $resultsPos =~ /\S/) {
	@FilesForMatrix = ($resultsNeg,$resultsPos);
		foreach $filereport (@FilesForMatrix) {
		chomp $filereport;
			if ($filereport =~ /(\S+\d+)(\/gsea_report_for_na_)([a-z]+)(_\d+)/) {
			$typeOfScore = $3;
			$baseFromReport = "$1";
			open REPORT, "<$filereport" or die "Can't open '$filereport'\n";
			$van = 0;
				while ($line = <REPORT>) {
				$van++;
					unless ($van == 1) {
					@arr = split ("\t", $line);
					$classid = @arr[0];
					$size = @arr[3];
					$es =  @arr[4];
					$nes =  @arr[5];
					$nompval = @arr[6];
					$fdrqval = @arr[7];
					
						### An rare exception was found where the filed corresponding to $nompval was empty in the output from gsea-3*
						### Since the FDR=1 I'm assigning a P-value=1 too
						if ($nompval =~ /^$/) {
							if ($fdrqval == 1) {
							$nompval = 1;
							}else{
							die "\n\nERROR!!! unexpected format in 'NOM p-val' = '$nompval' and 'FDR q-val' = '$fdrqval'\nin file:\n$filereport\nline:\n$line\n\n";
							}
						}
					
					$TypeOfScoreClass = "$typeOfScore\t$classid";
					$hashSummary{$ColumnName}{$TypeOfScoreClass}{nompval} = "_$nompval";
					$hashSummary{$ColumnName}{$TypeOfScoreClass}{fdrqval} = "_$fdrqval";
					$hashSummary{$ColumnName}{$TypeOfScoreClass}{real} = $size;
					$hashAllClassesTypes{$TypeOfScoreClass} = 1;
					
					### Here checking that enriched class specific files are available. If not, run the GSEA for this column again
						if ($nompval <= $hashParameters{cutoff_print_p} && $fdrqval <= $hashParameters{cutoff_print_q}) {
						$infileToIndex = "$baseFromReport/$classid.xls";
							if (-f $infileToIndex) {
							}else{
							$hashReRuns{$InFileName}{$classid} += 1;
								if ($hashReRuns{$InFileName}{$classid} == 1) {
								print "WARNING!!! couldn't get enriched class specific file:\n$infileToIndex\nP=$nompval; Q=$fdrqval\nRunning GSEA again with -plot_top_x = $PlotTopMany\n\n";
								&SendToGseaInBatchVerifyingOutfilesForGeneLevel($InFileName,$OutFileName,$PlotTopMany);
								}else{
								die "\n\nERROR!!! couldn't get enriched class specific file:\n$infileToIndex\nEven after running GSEA again with -plot_top_x = $PlotTopMany. Is the geneset below such thershold? If so increase -plot_top_x even more and try again\n\n";
								}
							}
						}
					}
				}
			close REPORT;
			}
		}
	}else{
	die "\n\nERROR!!! couldn't get reports for '$OutFileName'\n";
	}
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

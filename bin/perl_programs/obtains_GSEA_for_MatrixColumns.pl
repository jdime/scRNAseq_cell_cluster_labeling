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

$GseaExecutable =~ s/~\//\/$Users_home\/$DefaultUserName\//;

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
### Tables of P-value and FDR values obtained from GSEA for each column of -infile_matrix vs. each gene-class of -infile_classes,
### for 
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles    (path/name to the directory where outfiles will be saved)
###   -prefix_outfiles  (a string for the outfile names)
###   -infile_matrix    (path/name to the matrix of genes vs. conditions)
###   -infile_classes   (path/name to either the gene classes)
###   -nperm            (number of permutations to generate GSEA statistical analyses. GSEA developers recommend minimum '1000')
###   -use_universe     (indicates set of genes to use as universe from -infile_matrix and/or -infile_classes
###                      Either 'n1' to use all genes in -infile_matrix,
###                      Or 'i' to use the intersection of -infile_matrix and -infile_classes)
###   -cutoff_print_p   (once GSEA is computed, this cutoff selects Nom P-values that will be printed out in summary outfiles (e.g. 0.01))
###                      Note: also will generate outfiles without filtering neither P-values nor Q-values
###   -cutoff_print_q   (once GSEA is computed, this cutoff selects FDR Q-values that will be printed out in summary outfiles (e.g. 0.25))
###                      Note: also will generate outfiles without filtering neither P-values nor Q-values
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

### This will be used to get overall computing time
($year_overall1,$month_overall1,$day_overall1,$hour_overall1,$minute_overall1,$second_overall1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

&Readme;
&Parameters;

## By default GSEA plots only the top 20 sets (including detailed reports on 'core' gene sets, necessary to generate gene-level reports
## This may not be enough to gather all info needed for the full report on all enriched sets
## thus we'll first plot a few sets ($PlotTopFew) and if they are not enough then plot many ($PlotTopMany)
## In any case this shouldn't affect the class-level enrichment outfiles
$PlotTopFew  = 50;  
$PlotTopMany = 2000;

unless ($hashParameters{use_universe} =~ /^(i|n1)$/i) {
die "\n\nERROR!!! unexpected option '-use_universe $hashParameters{use_universe}'\n\n";
}

$useForLogOf0 = 0.0000000000001; ### To be use for -log(0) of p-values

############################
#### Step 2 -- Check that dependencies are available
############################

&CheckThatDependenciesAreAvailable();

############################
#### Step 3 -- Removing previous (presumably aborted) results
############################

$outdir = "$hashParameters{path_outfiles}/GSEA/$hashParameters{prefix_outfiles}";

if (-d $outdir) {
`rm -r $outdir`;
print "\n\nWARNING!!! removing preexisting outdir:\n'$outdir'\n";
}
system "mkdir -p $outdir";

############################
#### Step 4 -- Loading data and preparing infiles for GSEA
############################

### NOTE: load -infile_classes before -infile_matrix in case '-use_universe i' is used

$outfile_filtered_classes = "$hashParameters{path_outfiles}/$outfileWOpath_classes.Filtered.gmt";

LoadData::LoadClasses::LoadClasses($hashParameters{infile_classes},"ALL",1,1000000,"NA","NA","yes",$outfile_filtered_classes);
LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue::LoadMatrixReturnRnkEachColumnDepurateByMaxValue("$hashParameters{path_outfiles}/GSEA/$hashParameters{prefix_outfiles}",$hashParameters{infile_matrix},"NA","NA","NA",$hashParameters{use_universe},\%HashEachKeyClassesPassingCutoff);

############################
#### Step 5 -- Sending each column to GSEA
############################

print "Performing GSEA for each column\n";

### This will be used to get GSEA computing time
($year_gsea1,$month_gsea1,$day_gsea1,$hour_gsea1,$minute_gsea1,$second_gsea1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

foreach $c (1..$NumberOfColumnHeaders) {
$InputtedColumnName = $hashColumnsNames{$c};
$InFileName  = "$hashParameters{path_outfiles}/GSEA/$hashParameters{prefix_outfiles}/$InputtedColumnName.rnk";
$OutFileName = "$InputtedColumnName.vs.$outfileWOpath_classes.GSEA";

print "\nPerforming '$InputtedColumnName' of '$hashParameters{infile_matrix}' vs. $hashParameters{infile_classes} ($InFileName)'\n";
&SendToGseaInBatchVerifyingOutfilesForGeneLevel($InFileName,$OutFileName,$PlotTopFew);
}

### This will be used to get GSEA computing time
($year_gsea2,$month_gsea2,$day_gsea2,$hour_gsea2,$minute_gsea2,$second_gsea2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

############################
#### Step 6 -- Generate outfiles
############################

print "Generate outfiles\n";

$hashAllScoreSignsAndClasses{$ScoreSign}{$classid} = 1;
$hashAllScoreSigns{$ScoreSign} = 1;

### Open outfiles and print top-left corner string

foreach $ScoreSign (keys %hashAllScoreSigns) {
$outmatPAll     =  "$outdir/$hashParameters{prefix_outfiles}.GSEA.Pval.Unfiltered.$ScoreSign.mat.tsv";
$outmatPFil     =  "$outdir/$hashParameters{prefix_outfiles}.GSEA.Pval.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.$ScoreSign.mat.tsv";
$outmatFdrAll   =  "$outdir/$hashParameters{prefix_outfiles}.GSEA.Qval.Unfiltered.$ScoreSign.mat.tsv";
$outmatFdrFil   =  "$outdir/$hashParameters{prefix_outfiles}.GSEA.Qval.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}.$ScoreSign.mat.tsv";
$outmatPAllmLog =  "$outdir/$hashParameters{prefix_outfiles}.GSEA.Pval.Unfiltered.$ScoreSign.mLog.mat.tsv";

open $outmatPAll  ,   ">$outmatPAll"     or die "Can't open '$outmatPAll'\n";
open $outmatPFil  ,   ">$outmatPFil"     or die "Can't open '$outmatPFil'\n";
open $outmatFdrAll,   ">$outmatFdrAll"   or die "Can't open '$outmatFdrAll'\n";
open $outmatFdrFil,   ">$outmatFdrFil"   or die "Can't open '$outmatFdrFil'\n";
open $outmatPAllmLog, ">$outmatPAllmLog" or die "Can't open '$outmatPAllmLog'\n";

print $outmatPAll     "Pvalue.All";
print $outmatPFil     "Pvalue.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}";
print $outmatFdrAll   "FDRQvalue.All";
print $outmatFdrFil   "FDRQvalue.P$hashParameters{cutoff_print_p}.Q$hashParameters{cutoff_print_q}";
print $outmatPAllmLog "Pvalue.All.mLog";

	### Print column headers
	
	foreach $classid (sort keys %{$hashAllScoreSignsAndClasses{$ScoreSign}}) {
	print $outmatPAll     "\t$classid";
	print $outmatPFil     "\t$classid";
	print $outmatFdrAll   "\t$classid";
	print $outmatFdrFil   "\t$classid";
	print $outmatPAllmLog "\t$classid";
	}

	print $outmatPAll     "\n";
	print $outmatPFil     "\n";
	print $outmatFdrAll   "\n";
	print $outmatFdrFil   "\n";
	print $outmatPAllmLog "\n";
		
	### Print row headers and data

	foreach $c (1..$NumberOfColumnHeaders) {
		if ($hashColumnsNames{$c}) {
		$InputtedColumnName = $hashColumnsNames{$c};

		$ToPrintConcatenatedPvalAll     = "";
		$ToPrintConcatenatedQvalAll     = "";
		$ToPrintConcatenatedPvalFil     = "";
		$ToPrintConcatenatedQvalFil     = "";
		$ToPrintConcatenatedPvalAllmLog = "";
		
			foreach $classid (sort keys %{$hashAllScoreSignsAndClasses{$ScoreSign}}) {
			$PassPvalueFilter = 0;
			$PassQvalueFilter = 0;

				if ($hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{nompval}) {
				$nompval = $hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{nompval};
				$nompval =~ s/_//;
					if ($nompval <= $hashParameters{cutoff_print_p}) {
					$PassPvalueFilter = 1;
					}

					if ($nompval == 0) {
					$d = $useForLogOf0;
					}else{
					$d = $nompval;
					}
				$mLogNompval = log($d) * -1;
					
				}else{
				$nompval = "NA";
				$mLogNompval = "NA";
				}
	
				if ($hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{fdrqval}) {
				$fdrqval = $hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{fdrqval};
				$fdrqval =~ s/_//;
					if ($fdrqval <= $hashParameters{cutoff_print_q}) {
					$PassQvalueFilter = 1;
					}
				}else{
				$fdrqval = "NA";
				}

				if (($PassPvalueFilter == 1) && ($PassQvalueFilter == 1)) {
				$ToPrintConcatenatedPvalFil .= "\t$nompval";
				$ToPrintConcatenatedQvalFil .= "\t$fdrqval";
				}elsif ($nompval =~ /\d/) {
				$ToPrintConcatenatedPvalFil .= "\tNS";
				$ToPrintConcatenatedQvalFil .= "\tNS";
				}else{
				$ToPrintConcatenatedPvalFil .= "\tNA";
				$ToPrintConcatenatedQvalFil .= "\tNA";
				}
			$ToPrintConcatenatedPvalAll     .= "\t$nompval";
			$ToPrintConcatenatedQvalAll     .= "\t$fdrqval";
			$ToPrintConcatenatedPvalAllmLog .= "\t$mLogNompval";
			}
		$ToPrintConcatenatedPvalAll     =~ s/^\t//;
		$ToPrintConcatenatedQvalAll     =~ s/^\t//;
		$ToPrintConcatenatedPvalFil     =~ s/^\t//;
		$ToPrintConcatenatedQvalFil     =~ s/^\t//;
		$ToPrintConcatenatedPvalAllmLog =~ s/^\t//;
		
		print $outmatPAll     "$InputtedColumnName\t$ToPrintConcatenatedPvalAll\n";
		print $outmatFdrAll   "$InputtedColumnName\t$ToPrintConcatenatedQvalAll\n";
		print $outmatPFil     "$InputtedColumnName\t$ToPrintConcatenatedPvalFil\n";
		print $outmatFdrFil   "$InputtedColumnName\t$ToPrintConcatenatedQvalFil\n";
		print $outmatPAllmLog "$InputtedColumnName\t$ToPrintConcatenatedPvalAllmLog\n";
		}
	}
close $outmatPAll;
close $outmatPFil;
close $outmatFdrAll;
close $outmatFdrFil;
close $outmatPAllmLog;
}


$lookOrignalColumnNamesFile = "$outfileWOpath_matrix.OriginalColNames.txt";
if (-f $lookOrignalColumnNamesFile) {
`mv $lookOrignalColumnNamesFile $hashParameters{path_outfiles}`;
}

############################
#### Step 7 -- Generating log file with parameters used in the run
############################

&PrintParameters;

############################
#### Step 8 -- Remove temporary files
############################

print "Removing temporary files\n";

`rm $outfile_filtered_classes`;

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
$InputtedColumnName = $hashColumnsNames{$c};
$InFileName  = "$hashParameters{path_outfiles}/GSEA/$hashParameters{prefix_outfiles}/$InputtedColumnName.rnk";
$OutFileName = "$InputtedColumnName.vs.$outfileWOpath_classes.GSEA";
system "rm $InFileName";
system "rm -r $outdir/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/";
}

############################
#### Step 9 -- Report run time
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
#### Step 10 -- Finishing program
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

if ($hashParameters{infile_classes} =~ /-/) {
die "\n\nERROR!!! GSEA can't handle en dashes '-' in file names:\n'$hashParameters{infile_classes}'\n\n";
}
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_classes});
$outfileWOpath_classes = $outfileWOpath;

ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_matrix});
$outfileWOpath_matrix = $outfileWOpath;


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
$DirectoryFromPreviousRunLook = "$outdir/$OutFileName\.GseaPreranked\.[0-9]*[0-9]";

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
-gmx $outfile_filtered_classes \\
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
-out $outdir \\
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
$resultsNeg = "";
$resultsPos = "";
$resultsNeg = `ls $outdir/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/gsea_report_for_na_neg_[0-9]*[0-9]\.xls`;
$resultsPos = `ls $outdir/$OutFileName\.GseaPreranked\.[0-9]*[0-9]/gsea_report_for_na_pos_[0-9]*[0-9]\.xls`;

	if ($resultsNeg =~ /\S/ && $resultsPos =~ /\S/) {
	@FilesForMatrix = ($resultsNeg,$resultsPos);
		foreach $filereport (@FilesForMatrix) {
		chomp $filereport;
			if ($filereport =~ /(\S+\d+)(\/gsea_report_for_na_)([a-z]+)(_\d+)/) {
			$ScoreSign = $3;
			$baseFromReport = "$1";
			open REPORT, "<$filereport" or die "Can't open '$filereport'\n";
			$van = 0;
				while ($line = <REPORT>) {
				$van++;
					unless ($van == 1) {
					@arr = split ("\t", $line);
					$classid = @arr[0];
					$size    = @arr[3];
					$es      = @arr[4];
					$nes     = @arr[5];
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
					
					$hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{nompval} = "_$nompval";
					$hashSummary{$InputtedColumnName}{$ScoreSign}{$classid}{fdrqval} = "_$fdrqval";
					$hashAllScoreSignsAndClasses{$ScoreSign}{$classid} = 1;
					$hashAllScoreSigns{$ScoreSign} = 1;
					
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

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub CheckThatDependenciesAreAvailable {

### Checking that java v1.8 or v1.9 is available

print "Checking that java v1.8 or v1.9 is available\n";

	sub GetJavaVersion {
	$java_version = "";
	`java -version 2> $hashParameters{path_outfiles}/java_check.txt`;
	$java_version = `sed -n 1p $hashParameters{path_outfiles}/java_check.txt`;
	chomp $java_version;
	sleep 1;
	return $java_version;
	}


&GetJavaVersion();
if ($java_version =~ /(java|openjdk) version \"1\.[8|9]/) {
print "\tOK - found $java_version\n";
}else{
`module load java`;
sleep 5;
&GetJavaVersion();
	if ($java_version =~ /(java|openjdk) version \"1\.[8|9]/) {
	print "\tOK - found $java_version\n";
	}else{
	die "\n\nERROR!!! Couldn't find java (v1.8* or 1.9*) available\n"
	}

}

### Checking that GSEA is available

print "Checking that GSEA is available\n";

$GseaExecutable =~ s/~\//\/$Users_home\/$DefaultUserName\//;

### Checking if 'scratch' directory is available

$HomeUserDir = "\/$Users_home\/$DefaultUserName\/";
$ScratchUserDir = $HomeUserDir;
$ScratchUserDir =~ s/home/scratch/;

### Checking if 'scratch' contains a gsea-3.0.jar copy
### will use that one to write to cache files into scratch

if (-d $ScratchUserDir) {
$GseaExecutableInScratch = $GseaExecutable;
$GseaExecutableInScratch =~ s/home/scratch/;
	if (-f $GseaExecutableInScratch) {
	$GseaExecutable = $GseaExecutableInScratch;
	}
}

if (-f $GseaExecutable) {
print "\tOK - found $GseaExecutable\n";
}else{
die "\n\nERROR!!! Couldn't find gsea-3.0.jar executable at '$GseaExecutable'\n"
}

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

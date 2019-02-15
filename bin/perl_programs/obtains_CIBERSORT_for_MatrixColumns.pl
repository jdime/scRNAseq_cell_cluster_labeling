#!/usr/bin/perl

########################
### EXTERNAL DEPENDENCIES
########################
### (1)
### CIBERSORT.jar  -- (CIBERSORT standalone executable java *.jar file)
### It can be obtained from https://cibersort.stanford.edu/download.php
###
### In this script specify the path/name to the CIBERSORT.jar by variable $CibersortExecutable
### Make sure this file is at the same folder than the 'lib' folder provided by the Cibersort developers
###
### (2)
### transforms_GMT_to_Matrix_1_0.pl (only if using option '-classes_format gmt')
###
### (3)
### Java (tested with version 1.8.0_162)
###
### (4)
### R (tested with version 3.5.1) and libraries required by CIBERSORT ('e1071', 'parallel', 'preprocessCore', 'colorRamps', 'Rserve')
###
### (5)
### Perl modules -- Defined by 'use' command below
########################

############################
#### Step 1 -- Define parameters, default modules, infiles, and other dependencies
############################

use LoadParameters::Parameters;
use ReformatPerlEntities::ObtainOutfileWOpath;
use PathsDefinition::PathsToInputs;
use Date::Calc qw(Delta_DHMS);

### Put here the path/name of your CIBERSORT.jar executable
$CibersortExecutable = "~/PROGRAMS/CIBERSORT/SOFTWARE/CIBERSORT.jar";

$ThisProgramName = $0;
$ThisProgramName =~ s/\S+\///;

$CommentsForHelp = "
#####################################################################################
################### START INSTRUCTIONS TO RUN THIS PROGRAM ##########################
###
### Runs CIBERSORT for each column in -infile_matrix of genes (rows) vs. conditions (columns; e.g. samples, clusters, etc)
### and an -infile_classes with either cell-type marker gene classes or cell-type signatures (training gene expression profiles)
###
### Note: needs R library(Rserve) active
###       To activate it, in an R console type:
###       install.packages(\"Rserve\") ### only needed once
###       library(Rserve)
###       Rserve(args=\"--no-save\")
###
###       Check CIBERSORT's file CIBERSORT_documentation.txt for help
###
### -------------------------------------------INPUTS----------------------------------------------
###
### [1]
### a -infile_matrix in <tab> delimited format, like:
### GENE    cond_1  cond_2  cond_3...etc
### Gene_1  0.01    0.011   0.03
### Gene_2  -0.01   0.013   -0.13
### Gene_3  -0.02   0.014   -0.14
### ...etc
###
### [2a]
### a -infile_classes in <tab> delimited *.gmt format, like:
### CellType_1_id  CellType_1_name  Gene1  Gene2  Gene3
### CellType_2_id  CellType_2_name  Gene4  Gene5
### CellType_3_id  CellType_3_name  Gene5  Gene6
### ...etc
###
### OR
###
### [2b]
### a -infile_classes in <tab> delimited *.tsv format with gene expression profiles, like:
### GENE    CellType_1_id  CellType_2_id  CellType_3_id ... etc
### Gene_1  0.9            0.01           0.01
### Gene_2  0.8            0.02           0.02
### Gene_3  0.7            0.03           0.03
### Gene_4  0.02           0.8            0.04
### Gene_5  0.01           0.9            0.8
### Gene_6  0.01           0.01           0.9
### ...etc
###
### Notes:
###       a) if -infile_classes type 2a is provided, it will be converted to a type 2b-like binary matrix
###       b) -infile_classes type 2b is the 'signature' file used by CIBERSORT
###
### ----------------------------------------MAIN OUTPUTS-------------------------------------------
###
### [1]
### Tables of CIBERSORT enrichment scores for each column of -infile_matrix vs. each cell-type class or signature of -infile_classes
###
### [2]
### Tables of CIBERSORT p-values, Pearson correlation's and RMSE for each column of -infile_matrix
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles      (path/name to the directory where outfiles will be saved)
###   -prefix_outfiles    (a string for the outfile names)
###   -infile_matrix      (path/name to the matrix of genes vs. conditions)
###   -infile_classes     (path/name to either the cell-type marker gene classes or cell-type signatures)
###   -classes_type       (indicates if -infile_classes is either type 2a [use 'gmt'] or 2b [use 'profile'])
###   -nperm              (number of permutations to generate CIBERSORT statistical analysis. CIBERSORT developers recommend minimum '100')
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

### This will be used to get overall computing time
($year_overall1,$month_overall1,$day_overall1,$hour_overall1,$minute_overall1,$second_overall1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

&Readme;
&Parameters;

$CibersortExecutable =~ s/~\//\/$Users_home\/$DefaultUserName\//;

unless (-f $CibersortExecutable) {
die "\n\nERROR!!! Couldn't find CIBERSORT.jar executable at '$CibersortExecutable'\n"
}

############################
#### Step 1 -- Check that R(Rserve) is active
############################

print "Checking that R(Rserve) is active\n";

@Rserve = split ("\n", `ps ax | grep Rserve`);

$FoundRserveActive = 0;
foreach $rs (@Rserve) {
	if ($rs =~ /(\s*\d+\s+\?+\s+Ss\s+\d+:\d+)(\.\d+)*(\s+\S+)(Rserve)/) {
	$FoundRserveActive = 1;
	}
}
if ($FoundRserveActive == 1) {
print "OK - found Rserve active\n\n";
}else{
die "\n\nERROR!!! couldn't find R library(Rserve) active.
To activate it, in an R console type:
install.package(\"Rserve\") ### only needed once
library(Rserve)
Rserve(args=\"--no-save\")\n
Check CIBERSORT's file CIBERSORT_documentation.txt for help\n\n";
}
	

############################
#### Step 2 -- Getting original column headers to put them back in Cibersort's output
############################
@ColHeaders = split ("\t", `head -n 1 $hashParameters{infile_matrix}`);

$NumberOfColHeaders = -1;
foreach $colheader (@ColHeaders) {
chomp $colheader;
$NumberOfColHeaders++;
	unless ($NumberOfColHeaders == 0) {
	$hashColumnHeaders{$NumberOfColHeaders} = $colheader;
	}
}

############################
#### Step 3 -- Check format of -infile_classes
############################

if ($hashParameters{classes_type} =~ /^gmt$/i) {
&transforms_GMT_to_Matrix_1_0($hashParameters{infile_classes});
$ClassesFileForCibersort = "$hashParameters{path_outfiles}/$outfileWOpath_classes.1_0.tsv";
}else{
$ClassesFileForCibersort = "$hashParameters{infile_classes}";
}

############################
#### Step 4 -- Sending to CIBERSORT
############################

$outdir = "$hashParameters{path_outfiles}/CIBERSORT";
system "mkdir -p $outdir";

print "Performing CIBERSORT\n";

### This will be used to get CIBERSORT computing time
($year_cibersort1,$month_cibersort1,$day_cibersort1,$hour_cibersort1,$minute_cibersort1,$second_cibersort1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

system "date +%Y_%m_%d_%H_%M_%S";
$CommadToSendCibersort = "java -Xmx4g -Xms4g -jar $CibersortExecutable -M $hashParameters{infile_matrix} -B $ClassesFileForCibersort -n $hashParameters{nperm} > $outdir/$hashParameters{prefix_outfiles}.CIBERSORT.tmp";
print "\nRunning:\n'$CommadToSendCibersort'\n\n";
system "$CommadToSendCibersort";
system "date +%Y_%m_%d_%H_%M_%S";

### This will be used to get CIBERSORT computing time
($year_cibersort2,$month_cibersort2,$day_cibersort2,$hour_cibersort2,$minute_cibersort2,$second_cibersort2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

############################
#### Step 5 -- Reformat outfile headers
############################

open  OUTTMP, "<$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.tmp"                   or die "Couldn't open '$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.tmp'\n";
open  OUTES,  ">$outdir/$hashParameters{prefix_outfiles}.CIBERSORT_enrichment_scores.tsv" or die "Couldn't open '$outdir/$hashParameters{prefix_outfiles}.CIBERSORT_enrichment_scores.tsv'\n";
open  OUTPCR, ">$outdir/$hashParameters{prefix_outfiles}.CIBERSORT_PvalCorRmse.tsv"       or die "Couldn't open '$outdir/$hashParameters{prefix_outfiles}.CIBERSORT_PvalCorRmse.tsv'\n";

print "Reformat outfile headers and get final outfiles\n";

$firstOtherColumn = 0;
while ($line = <OUTTMP>) {
chomp $line;
$c = 0;
	if ($line =~ /^>/) {
	### omit
	}elsif ($line =~ /^(Column)/) {
	@ColHeaders = split ("\t", $line);
	print OUTES  "Input_Sample";
	print OUTPCR "Input_Sample";

		foreach $ch (@ColHeaders) {
		chomp $ch;
		$c++;
			if ($ch =~ /^P-value$/) {
			$firstOtherColumn = $c;
			}
			unless ($c == 1) {
				if ($firstOtherColumn > 0) {
				print OUTPCR "\t$ch";
				}else{
				print OUTES  "\t$ch";
				}
			}
		}
			
	print OUTES  "\n";
	print OUTPCR "\n";
		
	}else{
	@Data_es = split ("\t", $line);
	
		foreach $data (@Data_es) {
		$c++;
			if ($c == 1) {
			$ncolPlusOne = $data + 1;
				if ($hashColumnHeaders{$ncolPlusOne}) {
				$colheader = $hashColumnHeaders{$ncolPlusOne};
				print OUTPCR "$colheader";
				print OUTES  "$colheader";
				}else{
				die "\nERROR!!! couldn't find hashColumnHeaders{$ncolPlusOne}\n";
				}
			}elsif ($c >= $firstOtherColumn) {
			print OUTPCR "\t$data";
			}else{
			print OUTES  "\t$data";
			}
		}
			
	print OUTES  "\n";
	print OUTPCR "\n";
	}
}
close OUTES;

############################
#### Step 6 -- Generating log file with parameters used in the run
####           NOTE: this step must be run before getting gene-level networks
############################

&PrintParameters;

############################
#### Step 7 -- Remove temporary files
############################

print "Removing temporary files\n";

`rm $outdir/$hashParameters{prefix_outfiles}.CIBERSORT.tmp`;

############################
#### Step 8 -- Report run time
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

### cibersort runtime
($days_cibersort_diff, $hours_cibersort_diff, $minutes_cibersort_diff, $seconds_cibersort_diff) =
Delta_DHMS($year_cibersort1,$month_cibersort1,$day_cibersort1,$hour_cibersort1,$minute_cibersort1,$second_cibersort1,
	   $year_cibersort2,$month_cibersort2,$day_cibersort2,$hour_cibersort2,$minute_cibersort2,$second_cibersort2);
	   
$day_cibersort_seconds     = $days_cibersort_diff    * 86400;
$hour_cibersort_seconds    = $hours_cibersort_diff   * 3600;
$minute_cibersort_seconds  = $minutes_cibersort_diff * 60;
$all_cibersort_seconds     = $day_cibersort_seconds + $hour_cibersort_seconds + $minute_cibersort_seconds + $seconds_cibersort_diff;

open OUTRUNTIME, ">$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.Runtime" or die "Can't open '$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.Runtime'\n";
print OUTRUNTIME "cibersort\t$all_cibersort_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";
close OUTRUNTIME;

print "Took time:
cibersort\t$all_cibersort_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";

############################
#### Step 9 -- Finishing program
############################

print "\n  Done!!!\n\nCIBERSORT conducted for:\n'$hashParameters{infile_matrix}'\nvs.\n'$hashParameters{infile_classes}'\n\nCheck directory:\n$outdir\n\n";

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
'classes_type' => 1,
'infile_classes' => 1,
'nperm' => 1,
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
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_classes});
$outfileWOpath_classes = $outfileWOpath;

#### Ends -- Evaluate parameters
#######################

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub PrintParameters {
### Printing out parameters. Need to be concatenated at sub(Parameters)
open PARAMETERS, ">$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.Parameters" or die "Can't open '$outdir/$hashParameters{prefix_outfiles}.CIBERSORT.Parameters'\n";
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
sub transforms_GMT_to_Matrix_1_0 {

my($infile_gmt) = @_;

open INFILE_GMT, "<$infile_gmt" or die "Can't open -infile_gmt '$infile_gmt'\n\n";

while ($line = <INFILE_GMT>) {
chomp $line;

	unless ($line =~ /^#/) {
	@arr = split ("\t", $line);
	$c = 0;
		foreach $i (@arr) {
		$c++;
			if ($c == 1) {
			$classid = $i;
			}elsif ($c == 2) {
			### Do nothing
			}else{
			$hashDataForClasses{$classid}{$i} = 1;
			$hashAllGenes{$i} = 1;
			}
		}
	}
}
close INFILE_GMT;

open OUTFILE_BIN, ">$hashParameters{path_outfiles}/$outfileWOpath_classes.1_0.tsv" or die "Can't open '$hashParameters{path_outfiles}/$outfileWOpath_classes.1_0.tsv'\n\n";

print OUTFILE_BIN "GENES";

foreach $classid (sort keys %hashDataForClasses) {
print OUTFILE_BIN "\t$classid";
}
print OUTFILE_BIN "\n";

foreach $gene (sort keys %hashAllGenes) {
print OUTFILE_BIN "$gene";
	foreach $classid (sort keys %hashDataForClasses) {
		if ($hashDataForClasses{$classid}{$gene}) {
		print OUTFILE_BIN "\t1.0";
		}else{
		print OUTFILE_BIN "\t0.0";
		}
	}
print OUTFILE_BIN "\n";
}
close OUTFILE_BIN;

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#!/usr/bin/perl

########################
### See documentation below, in section 'START INSTRUCTIONS TO RUN THIS PROGRAM'
########################

########################
### EXTERNAL DEPENDENCIES
### R            -- functions ('fisher.test' and 'p.adjust')
### Perl modules -- Defined by 'use' command below
########################

use ReformatNumbers::NAvalues;
use LoadParameters::Parameters;
use LoadData::LoadMatrixNoStrain;
use LoadData::LoadClasses;
use ReformatPerlEntities::ObtainOutfileWOpath;
use PathsDefinition::PathsToInputs;
use Date::Calc qw(Delta_DHMS);

$ThisProgramName = $0;
$ThisProgramName =~ s/\S+\///;

$CommentsForHelp = "
#####################################################################################
################### START INSTRUCTIONS TO RUN THIS PROGRAM ##########################
###
### Runs Over Representation Analysis (ORA) for each column in -infile_matrix of genes (rows) vs. arrays (columns) e.g. samples, clusters, etc
### using an -infile_classes with gene sets (e.g. cell-type markers, Gene Ontology terms, etc)
### 
### To compute the ORA, uses a Fisher's exact test foreach column of -infile_matrix for each gene set of -infile_classes
###
### The ORA test is computed:
###  foreach column 'c' of -infile_matrix
###    foreach set/class 's' from -infile_classes
###      foreach cutoff sense 'e' (either -cutoff_neg or -cutoff_pos)
###
###      using the following contingency table which is passed to R function 'fisher.test'
###
###                             Genes passing 'e'  Genes not passing 'e'
###                               cutoff in 'c'       cutoff in 'c'
### Genes belonging to 's'             a                  b
### Genes not belonging 's'            c                  d
###
### From contingency table we get:
###       m = (a + b) = all genes belonging to a given class
###       n = (a + c) = all genes passing cutoff
###       o = (c + d) = all genes not belonging to a given class
###       p = (b + d) = all genes not passing cutoff
###       N = (a + b + c + d) = all genes in the universe
###
### To compute these fractions user can opt to use either the union or intersection of -infile_matrix and -infile_classes as the universe of genes
###
### Notes: 1) 'NA' values will be excluded
###
### -------------------------------------------INPUTS----------------------------------------------
###
### [1]
### -infile_matrix is the matrix of scores in <tab> delimited format:
### COMPLETE   Array_1   Array_2   Array_3
### Gene1      0.01      0.01      0.03
### Gene2      -0.01     NODATA    -0.15
### Gene3      NODATA    0.05      0.7
###
### [2] -infile_classes in *.gmt <tab> delimited format:
### CellType_1_id  CellType_1_name  Gene1  Gene2  Gene3
### CellType_2_id  CellType_2_name  Gene4  Gene5
### CellType_3_id  CellType_3_name  Gene5  Gene6
### ...etc
###
### ----------------------------------------MAIN OUTPUTS-------------------------------------------
###
### [1]
### Tables of p-values for each column header (from -infile_matrix) vs. each class (from -infile_classes)
### p-values are corrected for multiple hypothesis testing, and hence two files are provided (*Corrected* and *Uncorrected*)
###
### If requested by -cutoff_pos, overrepresentation analyzes are provided in *cutoff_pos* files
### If requested by -cutoff_neg, underrepresentation analyzes are provided in *cutoff_neg* files
###
### [2]
### -Log(p-value) is provided for the *ORA.Pval.Uncorrected.cutoff_pos.mat.txt
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles        (path/name to the directory where outfiles will be saved)
###   -prefix_outfiles      (string for outfiles' base name, or type 'default' to inherit from -infile_matrix and -infile_classes)
###
###   -infile_matrix        (path/name to the matrix with scores in <tab> delimited format)
###   -cutoff_neg           (cutoff to consider scores in -infile_matrix as negative-side 'white/drawn' balls for ORA, e.g. downregulated
###                          If your matrix only has positive scores set this parameter to 'NA')
###   -cutoff_pos           (cutoff to consider scores in -infile_matrix as positive-side 'white/drawn' balls for ORA, e.g. upregulated
###                          If your matrix only has negative scores set this parameter to 'NA')
###   -use_values_or_ranks  (indicates if cutoffs should be based on scores in -infile_matrix or their ranks
###                          For example, using '-use_values_or_ranks values -cutoff_pos 1.5' will consider values >= 1.5 as 'white/drawn' balls
###                          Whereas, using '-use_values_rank ranks -cutoff_pos 100' will consider the sorted top 100 hits as 'white/drawn' balls, independently of their value)
###
###   -infile_classes       (path/name to the list of class members in *gmt format)
###   -use_ora              (specifies if output will be filtered using either p-value [type 'p'] or p-value corrected [type 'pc'] from ORA test)
###   -cutoff_ora           (cutoff '<=' to consider a p-value from the ORA significant, e.g. '0.01' or '0.05')
###   -p_correction_test    (method for p-value correction, choose either 'bonferroni' (more stringent), 'BY' (medium stringent) or 'BH' (less stringent))
###   -use_universe         (indicates set of genes to use as universe from either union [type 'u'] or intersection [type 'i'] of -infile_matrix and -infile_classes)
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

### This will be used to get overall computing time
($year_overall1,$month_overall1,$day_overall1,$hour_overall1,$minute_overall1,$second_overall1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

&Readme;
&Parameters;

############################
#### Step 2 -- Loading data
############################


## Determining which cutoff types will be processed
if ($hashParameters{cutoff_neg} !~ /^NA$/i) { 
$ScoreSigns .= "\tcutoff_neg";
}
if ($hashParameters{cutoff_pos} !~ /^NA$/i) {
$ScoreSigns .= "\tcutoff_pos";
}
$ScoreSigns =~ s/^\t//;
@ScoreSigns = split ("\t", $ScoreSigns);

$outdir = "$hashParameters{path_outfiles}/ORA";
unless (-d $outdir) {
system "mkdir -p $outdir";
}

$useForLogOf0 = 0.0000000000001; ### To be use for -log(0) of p-values

### Hashing matrix
print "Reading Matrix\n";

LoadData::LoadMatrixNoStrain::LoadMatrixNoStrain($hashParameters{infile_matrix},"y","n");

print "Getting thersholds by ranks and values\n";

foreach $c (1..$NumberOfColumnHeaders) {
$NumberOfColumnHeaderToTrim = $c + 1;
$columnheader = $hashHeadersColumns{$c};

$command = "cut -f 1,$NumberOfColumnHeaderToTrim $hashParameters{infile_matrix} | sed -e 1d | sort -k 2 -g -r > $outdir/$outfileWOpath.$c.ForRank.tmp";
print "Column $c\n";
system $command;

## Get sorted list of row_id's and values
$NumberOfValues = `wc $outdir/$outfileWOpath.$c.ForRank.tmp`;
chomp $NumberOfValues;
$NumberOfValues =~ s/^\s+//;
$NumberOfValues =~ s/\s+.+//;
	
	unless ($hashParameters{cutoff_neg} =~ /^NA$/i) {
	$ThersholdRankNeg = $NumberOfValues - $hashParameters{cutoff_neg};
	}

	unless ($hashParameters{cutoff_pos} =~ /^NA$/i) {
	$ThersholdRankPos = $hashParameters{cutoff_pos};
	}

	open (RANKS, "<$outdir/$outfileWOpath.$c.ForRank.tmp") or die "Can't open '$outdir/$outfileWOpath.$c.ForRank.tmp'\n";
	$CountSortedValues = 0;
	while ($line = <RANKS>) {
	chomp $line;
	($rowheader,$value) = split ("\t", $line);

	$rowheader =~ tr/[a-z]/[A-Z]/;
	$value     =~ tr/[a-z]/[A-Z]/;
		
		unless ($hashNAvalues{$value}) {
		$CountSortedValues++;

		## here get list of rowheaders passing cutoffs by ranks and values
				
			unless ($hashParameters{cutoff_pos} =~ /^NA$/i) {
				if ($CountSortedValues <= $ThersholdRankPos) {
				$hashColumnRowPassingCutoffs{ranks}{cutoff_pos}{$columnheader}{$rowheader} = 1;
				#print "$columnheader\t$rowheader\t$CountSortedValues\n"; <STDIN>;
				}
				if ($value >= $hashParameters{cutoff_pos}) {
				$hashColumnRowPassingCutoffs{values}{cutoff_pos}{$columnheader}{$rowheader} = 1;
				}
			}
			unless ($hashParameters{cutoff_neg} =~ /^NA$/i) {
				if ($CountSortedValues >= $ThersholdRankNeg) {
				$hashColumnRowPassingCutoffs{ranks}{cutoff_neg}{$columnheader}{$rowheader} = 1;
				}
				if ($value <= $hashParameters{cutoff_neg}) {
				$hashColumnRowPassingCutoffs{values}{cutoff_neg}{$columnheader}{$rowheader} = 1;
				}
			}
		}
	}
	close RANKS;

system "rm $outdir/$outfileWOpath.$c.ForRank.tmp";
}


$outfile_filtered_classes = "$hashParameters{path_outfiles}/$outfileWOpath_classes.Filtered.gmt";

LoadData::LoadClasses::LoadClasses($hashParameters{infile_classes},"ALL",1,1000000,"NA","NA","yes",$outfile_filtered_classes);
system "rm $outfile_filtered_classes";

############################
#### Step 3 -- Get universe
############################

print "Get universe\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	foreach $r (1..$NumberOfRowHeaders) {
	$rowheader = $hashHeadersRows{$r};

		if ($hashData{$columnheader}{$rowheader}) {
		$score = $hashData{$columnheader}{$rowheader};
		$score =~ tr/[a-z]/[A-Z]/;
		$score =~ s/_//;
			
			unless ($hashNAvalues{$score}) {
			$hashUniverse{n1}{$columnheader}{$rowheader} = 1;
			$hashUniverse{u}{$columnheader}{$rowheader} = 1;
			}
		}
	}
}

foreach $rowheader (keys %HashEachKeyClassesPassingCutoff) {
	foreach $c (1..$NumberOfColumnHeaders) {
	$columnheader = $hashHeadersColumns{$c};

	$hashUniverse{n2}{$columnheader}{$rowheader} = 1;
	$hashUniverse{u}{$columnheader}{$rowheader} = 1;
		if ($hashUniverse{n1}{$columnheader}{$rowheader}) {
		$hashUniverse{i}{$columnheader}{$rowheader} = 1;
		}
	}
}

############################
#### Step 4 -- Compiling data to obtain ORA
############################

### This will be used to get HYPERG computing time
($year_hyperg1,$month_hyperg1,$day_hyperg1,$hour_hyperg1,$minute_hyperg1,$second_hyperg1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

print "Compiling data to obtain ORA\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	### Here restricting to selected universe
	foreach $rowheader (keys %{$hashUniverse{$hashParameters{use_universe}}{$columnheader}}) {
	
		if ($HashEachKeyClasses{$rowheader}) {
		@AllClassesThisGene = split (",", $HashEachKeyClasses{$rowheader});
		}

		### Indexing 'N'
		$hashDataForHyperG{$columnheader}{N} += 1;
				
		unless ($hashParameters{cutoff_pos} =~ /^na$/i) {

			### Indexing 'n'

			if ($hashColumnRowPassingCutoffs{$hashParameters{use_values_or_ranks}}{cutoff_pos}{$columnheader}{$rowheader}) {
			$hashDataForHyperG{$columnheader}{cutoff_pos}{n} += 1;
			
				### Indexing 'a'
				foreach $classid (@AllClassesThisGene) {

					if ($classid =~ /\S/) {
						if ($hashClassesRetained{$classid}) {
						$hashCountGenesInEachColumnEachClass{$columnheader}{$classid} += 1;
						$hashAllClassesInThisMatrix{$classid} = 1;
						$hashDataForHyperG{$columnheader}{cutoff_pos}{$classid}{a} += 1;
						$hashDataForGeneLevelNetwork{$columnheader}{cutoff_pos}{$classid} .= ",$rowheader";
						}
					}
				}
			}else{
			
			### Indexing 'p'
			
			$hashDataForHyperG{$columnheader}{cutoff_pos}{p} += 1;
			
				### Indexing 'b'
				foreach $classid (@AllClassesThisGene) {
					if ($classid =~ /\S/) {
						if ($hashClassesRetained{$classid}) {
						$hashCountGenesInEachColumnEachClass{$columnheader}{$classid} += 1;
						$hashAllClassesInThisMatrix{$classid} = 1;
						$hashDataForHyperG{$columnheader}{cutoff_pos}{$classid}{b} += 1;
						}
					}
				}
			}
		}
		
		unless ($hashParameters{cutoff_neg} =~ /^na$/i) {

			### Indexing 'n'

			if ($hashColumnRowPassingCutoffs{$hashParameters{use_values_or_ranks}}{cutoff_neg}{$columnheader}{$rowheader}) {
			$hashDataForHyperG{$columnheader}{cutoff_neg}{n} += 1;
			
				### Indexing 'a'
				foreach $classid (@AllClassesThisGene) {
					if ($classid =~ /\S/) {
						if ($hashClassesRetained{$classid}) {
						$hashCountGenesInEachColumnEachClass{$columnheader}{$classid} += 1;
						$hashAllClassesInThisMatrix{$classid} = 1;
						$hashDataForHyperG{$columnheader}{cutoff_neg}{$classid}{a} += 1;
						$hashDataForGeneLevelNetwork{$columnheader}{cutoff_neg}{$classid} .= ",$rowheader";
						}
					}
				}
			}else{
			
			### Indexing 'p'
			
			$hashDataForHyperG{$columnheader}{cutoff_neg}{p} += 1;
			
				### Indexing 'b'
				foreach $classid (@AllClassesThisGene) {
					if ($classid =~ /\S/) {
						if ($hashClassesRetained{$classid}) {
						$hashCountGenesInEachColumnEachClass{$columnheader}{$classid} += 1;
						$hashAllClassesInThisMatrix{$classid} = 1;
						$hashDataForHyperG{$columnheader}{cutoff_neg}{$classid}{b} += 1;
						}
					}
				}
			}
		}
	}
}

############################
#### Step 5 -- Preparing infiles to obtain ORA and sending them to R
############################
print "Preparing infiles to obtain ORA and sending them to R\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	### Obtaining 'N'
	if ($hashDataForHyperG{$columnheader}{N}) {
	$N = $hashDataForHyperG{$columnheader}{N};
	}else{
	die "ERROR!!! couldn't obtain hashDataForHyperG{$columnheader}{N}\n";
	}

	foreach $ScoreSign (@ScoreSigns) {

		unless ($hashParameters{$ScoreSign} =~ /^na$/i) {
		
		###Obtaining 'n' and 'p'
			if ($hashDataForHyperG{$columnheader}{$ScoreSign}{n}) {
			$n = $hashDataForHyperG{$columnheader}{$ScoreSign}{n};
			}else{
			$n = 0;
			print "WARNING!!! in column '$columnheader' none value passed $ScoreSign '$hashParameters{$ScoreSign}'\n";
			}

			if ($hashDataForHyperG{$columnheader}{$ScoreSign}{p}) {
			$p = $hashDataForHyperG{$columnheader}{$ScoreSign}{p};
			}else{
			$p = 0;
			print "WARNING!!! in column '$columnheader' all values passed $ScoreSign '$hashParameters{$ScoreSign}'\n";
			}
			
			$classesSentToR = 0;
			
			unless (($n == 0) or ($p == 0)) {

			#############################
			## Printing data matrix for R
			
			print "Printing data matrix for R\n\n";
		
			open MATFORR, ">$outdir/$outfileWOpath.$columnheader.$ScoreSign.matforR.mat" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$ScoreSign.matforR.mat'\n";
			print MATFORR "Column_header---Cutoff_type---Class_id\ta\tb\tc\td\n";

				foreach $classid (sort keys %hashAllClassesInThisMatrix) {
				$classesSentToR++;
				$hashClassesSentToRPositionName{$classesSentToR} = $classid;
				$hashClassesSentToRName{$classid} = 1;
						
					###Obtaining 'a' and 'b'
					if ($hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{a}) {
					$a = $hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{a};
					}else{
					$a = 0;
					}
	
					if ($hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{b}) {
					$b = $hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{b};
					}else{
					$b = 0;
					}
							
						
				###Obtaining 'c', 'd', 'm' and 'o'
						
				$c = $n - $a;
				$d = $p - $b;
				$m = $a + $b;
				$o = $c + $d;
					
				$ratio = ($a / $n) / ($m / $N);
						
				$hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{c} = "_$c";
				$hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{d} = "_$d";
				$hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{m} = "_$m";
				$hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{o} = "_$o";
				$hashDataForHyperG{$columnheader}{$ScoreSign}{$classid}{ratio} = "_$ratio";
					
					if ($hashMaxPopulationEachClassEachScoreSign{$classid}{$ScoreSign}) {
						if ($a > $hashMaxPopulationEachClassEachScoreSign{$classid}{$ScoreSign}) {
						$hashMaxPopulationEachClassEachScoreSign{$classid}{$ScoreSign} = $a;
						}
					}else{
					$hashMaxPopulationEachClassEachScoreSign{$classid}{$ScoreSign} = $a;
					}

				print MATFORR "$columnheader---$ScoreSign---$classid\t$a\t$b\t$c\t$d\n";
				}
				
			close MATFORR;
			
			###########################
			## Printing commands for R

			print "Printing/send commands for/to R\n\n";
			
			open INSFORR, ">$outdir/$outfileWOpath.$columnheader.$ScoreSign.InsFor.R" or die "Can't open '$outfileWOpath.$columnheader.$ScoreSign.InsFor.R'\n";
			print INSFORR "##Computing fisher.test for classes ocurring in $columnheader---$ScoreSign\n";
			print INSFORR "mat<-read.table(\"$outdir/$outfileWOpath.$columnheader.$ScoreSign.matforR.mat\",header=TRUE,row.names=1)
                                       dim(mat)
                                       NumberOfClassesToTest<-nrow(mat)
                                       NumberOfClassesToTest
                                       TableForOutfile.df <-data.frame(p.val.uncorrected=numeric(NumberOfClassesToTest),stringsAsFactors=FALSE)
                                       for (i in 1: NumberOfClassesToTest) {
                                        ContingencyTable<-matrix(data=c(mat\[i,\"a\"\],mat\[i,\"b\"\],mat\[i,\"c\"\],mat\[i,\"d\"\]),nrow=2,ncol=2)
                                        p.value.uncorrected<-fisher.test (ContingencyTable, alternative = \"greater\")\$p.value
                                        TableForOutfile.df\$p.val.uncorrected[i]<-(p.value.uncorrected)
                                       }
				       p.val.uncorrected<-TableForOutfile.df\$p.val.uncorrected
				       p.val.corrected<-p.adjust(p.val.uncorrected,method=\"$hashParameters{p_correction_test}\",n=NumberOfClassesToTest)
                                       write.table(TableForOutfile.df\$p.val.uncorrected,\"$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalUncorrected\",sep=\"\\t\",quote=FALSE,col.names=NA,row.names=TRUE)
                                       write.table(p.val.corrected,\"$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalCorrected\",sep=\"\\t\",quote=FALSE,col.names=NA,row.names=TRUE)
				       q()\n";
			close INSFORR;
				       
			###########################
			## Sending to R
			system "R --no-save < $outdir/$outfileWOpath.$columnheader.$ScoreSign.InsFor.R";

			###########################
			## Obtaining P-values and reformating matrices
			
			print "Obtaining P-values and reformating matrices\n\n";
			
			open UNCORRECTED, "<$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalUncorrected" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalUncorrected'\n";
			$linesHyperUnCorrected = 0;
				while ($line = <UNCORRECTED>) {
					if ($line =~ /^(\S+)(\t)(\S+)/) {
					$Puncorrected = $3;
					$linesHyperUnCorrected++;
					$classid = $hashClassesSentToRPositionName{$linesHyperUnCorrected};
					$hashFinalData{$columnheader}{$ScoreSign}{$classid}{puncorrected} = "_$Puncorrected";
					}
				}
			close UNCORRECTED;
				unless ($linesHyperUnCorrected == $classesSentToR) {
				die "ERROR!!! $classesSentToR classes were sent to R for ORA, but $linesHyperUnCorrected uncorrected P-values were retrieved\n";
				}
			
			open CORRECTED, "<$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalCorrected" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalCorrected'\n";
			$linesHyperCorrected = 0;
				while ($line = <CORRECTED>) {
					if ($line =~ /^(\S+)(\t)(\S+)/) {
					$Pcorrected = $3;
					$linesHyperCorrected++;
					$classid = $hashClassesSentToRPositionName{$linesHyperCorrected};
					$hashFinalData{$columnheader}{$ScoreSign}{$classid}{pcorrected} = "_$Pcorrected";
					}
				}
			close CORRECTED;
				unless ($linesHyperCorrected == $classesSentToR) {
				die "ERROR!!! $classesSentToR classes were sent to R for ORA, but $linesHyperCorrected corrected P-values were retrieved\n";
				}
				
			system "rm $outdir/$outfileWOpath.$columnheader.$ScoreSign.matforR.mat $outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalUncorrected $outdir/$outfileWOpath.$columnheader.$ScoreSign.PvalCorrected $outdir/$outfileWOpath.$columnheader.$ScoreSign.InsFor.R";
			
			}
		}
	}
}

### This will be used to get HYPERG computing time
($year_hyperg2,$month_hyperg2,$day_hyperg2,$hour_hyperg2,$minute_hyperg2,$second_hyperg2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);


############################
#### Step 6 -- Generating outfiles
############################
print "Generating final report\n\n";

foreach $ScoreSign (@ScoreSigns) {
	unless ($hashParameters{$ScoreSign} =~ /^na$/i) {
	$outmatPUnAll     =  "$outdir/$hashParameters{prefix_outfiles}.ORA.Pval.Uncorrected.$ScoreSign.mat.tsv";
	$outmatPUnFil     =  "$outdir/$hashParameters{prefix_outfiles}.ORA.Pval.Uncorrected.P$hashParameters{cutoff_ora}.$ScoreSign.mat.tsv";
	$outmatPCoAll     =  "$outdir/$hashParameters{prefix_outfiles}.ORA.Pval.Corrected.$ScoreSign.mat.tsv";
	$outmatPCoFil     =  "$outdir/$hashParameters{prefix_outfiles}.ORA.Pval.Corrected.P$hashParameters{cutoff_ora}.$ScoreSign.mat.tsv";
	$outmatPUnAllmLog =  "$outdir/$hashParameters{prefix_outfiles}.ORA.Pval.Uncorrected.$ScoreSign.mLog.mat.tsv";

	open $outmatPUnAll,     ">$outmatPUnAll"     or die "Can't open '$outmatPUnAll'\n";
	open $outmatPUnFil,     ">$outmatPUnFil"     or die "Can't open '$outmatPUnFil'\n";
	open $outmatPCoAll,     ">$outmatPCoAll"     or die "Can't open '$outmatPCoAll'\n";
	open $outmatPCoFil,     ">$outmatPCoFil"     or die "Can't open '$outmatPCoFil'\n";
	open $outmatPUnAllmLog, ">$outmatPUnAllmLog" or die "Can't open '$outmatPUnAllmLog'\n";

	print $outmatPUnAll     "ORA.Pval.Uncorrected.$ScoreSign";
	print $outmatPUnFil     "ORA.Pval.Uncorrected.P$hashParameters{cutoff_ora}.$ScoreSign";
	print $outmatPCoAll     "ORA.Pval.Corrected.$ScoreSign";
	print $outmatPCoFil     "ORA.Pval.Corrected.P$hashParameters{cutoff_ora}.$ScoreSign";
	print $outmatPUnAllmLog "ORA.Pval.Uncorrected.$ScoreSign.mLog";

		### Printing column headers
		
		foreach $classid (sort keys %hashClassesSentToRName) {
		print $outmatPUnAll     "\t$classid";
		print $outmatPUnFil     "\t$classid";
		print $outmatPCoAll     "\t$classid";
		print $outmatPCoFil     "\t$classid";
		print $outmatPUnAllmLog "\t$classid";
		}

		print $outmatPUnAll     "\n";
		print $outmatPUnFil     "\n";
		print $outmatPCoAll     "\n";
		print $outmatPCoFil     "\n";
		print $outmatPUnAllmLog "\n";


		### Print row headers and data
		
		foreach $c (1..$NumberOfColumnHeaders) {
		$columnheader = $hashHeadersColumns{$c};

		print $outmatPUnAll     "$columnheader";
		print $outmatPUnFil     "$columnheader";
		print $outmatPCoAll     "$columnheader";
		print $outmatPCoFil     "$columnheader";
		print $outmatPUnAllmLog "$columnheader";
			
			foreach $classid (sort keys %hashClassesSentToRName) {
	
				if ($hashFinalData{$columnheader}{$ScoreSign}{$classid}{puncorrected}) {
				$puncorrected = $hashFinalData{$columnheader}{$ScoreSign}{$classid}{puncorrected};
				$puncorrected =~ s/_//;
				}else{
				$puncorrected = "NA";
				}
				if ($hashFinalData{$columnheader}{$ScoreSign}{$classid}{pcorrected}) {
				$pcorrected = $hashFinalData{$columnheader}{$ScoreSign}{$classid}{pcorrected};
				$pcorrected =~ s/_//;
				}else{
				$pcorrected = "NA";
				}
				
				### Determining if passes the P-value cutoff
				$passesPvalCutoff = 0;
				if ($puncorrected =~ /^na$/i) {
				print $outmatPUnAll     "\t$puncorrected";
				print $outmatPUnFil     "\t$puncorrected";
				print $outmatPCoAll     "\t$pcorrected";
				print $outmatPCoFil     "\t$pcorrected";
				print $outmatPUnAllmLog "\t$puncorrected";
				
				}else{
					if ($hashParameters{use_ora} =~ /^pc$/) {
						if ($pcorrected <= $hashParameters{cutoff_ora}) {
						$passesPvalCutoff = 1;
						}
					}elsif ($hashParameters{use_ora} =~ /^p$/) {
						if ($puncorrected <= $hashParameters{cutoff_ora}) {
						$passesPvalCutoff = 1;
						}
					}
					
					if ($puncorrected == 0) {
					$d = $useForLogOf0;
					}else{
					$d = $puncorrected;
					}
					$mLogPuncorrected = log($d) * -1;
					
					if ($passesPvalCutoff == 1) {
					print $outmatPUnAll     "\t$puncorrected";
					print $outmatPUnFil     "\t$puncorrected";
					print $outmatPCoAll     "\t$pcorrected";
					print $outmatPCoFil     "\t$pcorrected";
					print $outmatPUnAllmLog "\t$mLogPuncorrected";
						
					}else{
					print $outmatPUnAll     "\t$puncorrected";
					print $outmatPUnFil     "\tNS";
					print $outmatPCoAll     "\t$pcorrected";
					print $outmatPCoFil     "\tNS";
					print $outmatPUnAllmLog "\t$mLogPuncorrected";
					}
				}
			}
		print $outmatPUnAll     "\n";
		print $outmatPUnFil     "\n";
		print $outmatPCoAll     "\n";
		print $outmatPCoFil     "\n";
		print $outmatPUnAllmLog "\n";
		}
	close $outmatPUnAll;
	close $outmatPUnFil;
	close $outmatPCoAll;
	close $outmatPCoFil;
	close $outmatPUnAllmLog;
	}
}

&Readme;
&PrintParameters;

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

### hyperg runtime
($days_hyperg_diff, $hours_hyperg_diff, $minutes_hyperg_diff, $seconds_hyperg_diff) =
Delta_DHMS($year_hyperg1,$month_hyperg1,$day_hyperg1,$hour_hyperg1,$minute_hyperg1,$second_hyperg1,
	   $year_hyperg2,$month_hyperg2,$day_hyperg2,$hour_hyperg2,$minute_hyperg2,$second_hyperg2);
	   
$day_hyperg_seconds     = $days_hyperg_diff    * 86400;
$hour_hyperg_seconds    = $hours_hyperg_diff   * 3600;
$minute_hyperg_seconds  = $minutes_hyperg_diff * 60;
$all_hyperg_seconds     = $day_hyperg_seconds + $hour_hyperg_seconds + $minute_hyperg_seconds + $seconds_hyperg_diff;

open OUTRUNTIME, ">$outdir/$outfileWOpath.ORA.Runtime" or die "Can't open '$outdir/$outfileWOpath.ORA.Runtime'\n";
print OUTRUNTIME "ORA\t$all_hyperg_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";
close OUTRUNTIME;

print "Took time:
ORA\t$all_hyperg_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";

print "\n  Done!!!\n\nCheck outfiles '$outdir/$outfileWOpath.ORA.Pval*'\n";

exit;
#############################################################
###################### END OF PROGRAM  ######################
#############################################################


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
'prefix_outfiles' => 1,
'infile_matrix' => 1,
'cutoff_neg' => 1,
'cutoff_pos' => 1,
'use_values_or_ranks' => 1,
'infile_classes' => 1,
'use_ora' => 1,
'cutoff_ora' => 1,
'p_correction_test' => 1,
'use_universe' => 1,
);

#######################
#### Starts -- Evaluate parameters

&LoadParameters::Parameters::MainSubParameters(\%hashParametersTolookFor,\@arrayInputtedOneLineCommands);
$Parameters .= "$MoreParameters";

ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_classes});
$outfileWOpath_classes = $outfileWOpath;

## Defining prefix string for OUTFILE

if ($hashParameters{prefix_outfiles} =~ /^default$/) {

ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_matrix});
$outfileWOpathA = $outfileWOpath;

ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile_classes});
$outfileWOpathB = $outfileWOpath;

$outfileWOpath = "$outfileWOpathA.$hashParameters{use_values_or_ranks}.CutoffNeg$hashParameters{cutoff_neg}.CutoffPos$hashParameters{cutoff_pos}.vs.$outfileWOpathB";

}else{
$outfileWOpath = $hashParameters{prefix_outfiles};
}


#### Ends -- Evaluate parameters
#######################

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub PrintParameters {

### Printing out parameters. Need to be concatenated at sub(Parameters)
open PARAMETERS, ">$outdir/$outfileWOpath.ORA.Parameters" or die "Can't open '$outdir/$outfileWOpath.ORA.Parameters'\n";
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

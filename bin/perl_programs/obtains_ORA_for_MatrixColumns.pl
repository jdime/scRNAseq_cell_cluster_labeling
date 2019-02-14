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
### Gene2     -0.01      NODATA   -0.15
### Gene3     NODATA     0.05      0.7
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
### Tables of column headers (from -infile_matrix) vs. significantly enriched classes (given a provided P-value cutoff)
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles        (path/name to the directory where outfiles will be saved)
###   -prefix_outfiles      (string for outfiles' base name, or type 'default' to inherit from -infile_matrix and -infile_classes)
###
###   -infile_matrix        (path/name to the matrix with scores in <tab> delimited format)
###   -numb_cols            (indicates columns from -infile_matrix to be processed e.g:
###                          to process initial column type '1'
###                          to process first two columns type '1,2'
###                          to process columns 1 to 4 and 7 type '1-4,7'
###                          to process all columns type 'ALL'
###   -cutoff_neg           (cutoff to consider scores in -infile_matrix as negative-side 'white/drawn' balls for ORA, e.g. downregulated
###                          If your matrix only has positive scores set this parameter to 'NA')
###   -cutoff_pos           (cutoff to consider scores in -infile_matrix as positive-side 'white/drawn' balls for ORA, e.g. upregulated
###                          If your matrix only has negative scores set this parameter to 'NA')
###   -use_values_or_ranks  (indicates if cutoffs should be based on scores in -infile_matrix or their ranks
###                          For example, using '-use_values_or_ranks values -cutoff_pos 1.5' will consider values >= 1.5 as 'white/drawn' balls
###                          Whereas, using '-use_values_rank ranks -cutoff_pos 100' will consider the sorted top 100 hits as 'white/drawn' balls, independently of their value)
###
###   -infile_classes       (path/name to the list of class members in *gmt format)
###   -restrict_classes     (one class identifier per-line. Or type 'ALL' to include all classes in -infile_classes. In any case population cutoffs -set_min and -set_max will be applied)
###   -set_min              (minimum number of genes in a class to consider such class in ORA analysis)
###   -set_max              (maximum number of genes in a class to consider such class in ORA analysis)
###
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

## Determining which cutoff types will be processed
if ($hashParameters{cutoff_neg} !~ /^NA$/i) { 
$cutoffTypes .= "\tcutoff_neg";
}
if ($hashParameters{cutoff_pos} !~ /^NA$/i) {
$cutoffTypes .= "\tcutoff_pos";
}
$cutoffTypes =~ s/^\t//;
@CutoffTypes = split ("\t", $cutoffTypes);

$outdir = "$hashParameters{path_outfiles}/ORA";
unless (-d $outdir) {
system "mkdir -p $outdir";
}

############################
#### Step 1 -- Loading data
############################
### Hashing matrix
print "Reading Matrix\n";

LoadData::LoadMatrixNoStrain::LoadMatrixNoStrain($hashParameters{infile_matrix},"y","n");

print "Getting thersholds by ranks and values\n";

foreach $c (1..$NumberOfColumnHeaders) {
$NumberOfColumnHeaderToTrim = $c + 1;
$columnheader = $hashHeadersColumns{$c};

	## here restricting to requested columns
	if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {
	$command = "cut -f 1,$NumberOfColumnHeaderToTrim $hashParameters{infile_matrix} | sed -e 1d | sort -k 2 -g -r > $outdir/$outfileWOpath.$c.ForRank.tmp";
	print "Column $c\n";
	system $command;

	## here get sorted list of row_id's and values
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
}

### Note although here I'm setting 'set_min' and 'set_max' is just to save computing time to avoid indexing very small/large classes,
### the actual filter by population will be applied below restricting to genes occurring in each column of this matrix, foreach class
LoadData::LoadClasses::LoadClasses($hashParameters{infile_classes},$hashParameters{restrict_classes},$hashParameters{set_min},$hashParameters{set_max},"NA","NA","yes");
system "rm $hashParameters{infile_classes}.Filtered.gmt";

############################
#### Step 2 -- Get universe
############################

print "Get universe\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	## here restricting to requested columns
	if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {

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
}

foreach $rowheader (keys %HashEachKeyClassesPassingCutoff) {
	foreach $c (1..$NumberOfColumnHeaders) {
	$columnheader = $hashHeadersColumns{$c};

		## here restricting to requested columns
		if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {
		$hashUniverse{n2}{$columnheader}{$rowheader} = 1;
		$hashUniverse{u}{$columnheader}{$rowheader} = 1;
			if ($hashUniverse{n1}{$columnheader}{$rowheader}) {
			$hashUniverse{i}{$columnheader}{$rowheader} = 1;
			}
		}
	}
}

############################
#### Step 3 -- Compiling data to obtain ORA
############################

### This will be used to get HYPERG computing time
($year_hyperg1,$month_hyperg1,$day_hyperg1,$hour_hyperg1,$minute_hyperg1,$second_hyperg1) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);

print "Compiling data to obtain ORA\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	## here restricting to requested columns
	if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {

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
}

############################
#### Step 4 -- Preparing infiles to obtain ORA and sending them to R
############################
print "Preparing infiles to obtain ORA and sending them to R\n\n";

foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};

	## here restricting to requested columns
	if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {

	### Obtaining 'N'
		if ($hashDataForHyperG{$columnheader}{N}) {
		$N = $hashDataForHyperG{$columnheader}{N};
		}else{
		die "ERROR!!! couldn't obtain hashDataForHyperG{$columnheader}{N}\n";
		}

		foreach $cutofftype (@CutoffTypes) {

			unless ($hashParameters{$cutofftype} =~ /^na$/i) {
			
			###Obtaining 'n' and 'p'
				if ($hashDataForHyperG{$columnheader}{$cutofftype}{n}) {
				$n = $hashDataForHyperG{$columnheader}{$cutofftype}{n};
				}else{
				$n = 0;
				print "WARNING!!! in column '$columnheader' none value passed $cutofftype '$hashParameters{$cutofftype}'\n";
				}

				if ($hashDataForHyperG{$columnheader}{$cutofftype}{p}) {
				$p = $hashDataForHyperG{$columnheader}{$cutofftype}{p};
				}else{
				$p = 0;
				print "WARNING!!! in column '$columnheader' all values passed $cutofftype '$hashParameters{$cutofftype}'\n";
				}
				
				$classesSentToR = 0;
				
				unless (($n == 0) or ($p == 0)) {

				#############################
				## Printing data matrix for R
				
				print "Printing data matrix for R\n\n";
			
				open MATFORR, ">$outdir/$outfileWOpath.$columnheader.$cutofftype.matforR.mat" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$cutofftype.matforR.mat'\n";
				print MATFORR "Column_header---Cutoff_type---Class_id\ta\tb\tc\td\n";

					foreach $classid (sort keys %hashAllClassesInThisMatrix) {
						if ($hashCountGenesInEachColumnEachClass{$columnheader}{$classid}) {
						$populationInThisColumn = $hashCountGenesInEachColumnEachClass{$columnheader}{$classid};
							
						### Here filtering classes by population in each column
							if (($populationInThisColumn >= $hashParameters{set_min}) && ($populationInThisColumn <= $hashParameters{set_max})) {
							
							$classesSentToR++;
							$hashClassesSentToRPositionName{$classesSentToR} = $classid;
							$hashClassesSentToRName{$classid} = 1;
							
							### Indexing maximum population of each class in this matrix (for print out below)
								if ($hashMaxPopulationEachClass{$classid}) {
								$populationInThisColumn > $hashMaxPopulationEachClass{$classid};
								$hashMaxPopulationEachClass{$classid} = $populationInThisColumn;
								}else{
								$hashMaxPopulationEachClass{$classid} = $populationInThisColumn;
								}
							
							###Obtaining 'a' and 'b'
								if ($hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{a}) {
								$a = $hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{a};
								}else{
								$a = 0;
								}
		
								if ($hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{b}) {
								$b = $hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{b};
								}else{
								$b = 0;
								}
								
							
							###Obtaining 'c', 'd', 'm' and 'o'
							
							$c = $n - $a;
							$d = $p - $b;
							$m = $a + $b;
							$o = $c + $d;
							
							$ratio = ($a / $n) / ($m / $N);
							
							$hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{c} = "_$c";
							$hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{d} = "_$d";
							$hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{m} = "_$m";
							$hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{o} = "_$o";
							$hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{ratio} = "_$ratio";
							
								if ($hashMaxPopulationEachClassEachCutofftype{$classid}{$cutofftype}) {
									if ($a > $hashMaxPopulationEachClassEachCutofftype{$classid}{$cutofftype}) {
									$hashMaxPopulationEachClassEachCutofftype{$classid}{$cutofftype} = $a;
									}
								}else{
								$hashMaxPopulationEachClassEachCutofftype{$classid}{$cutofftype} = $a;
								}
	
							print MATFORR "$columnheader---$cutofftype---$classid\t$a\t$b\t$c\t$d\n";
							
							}
						}
					}
				close MATFORR;
				
				###########################
				## Printing commands for R

				print "Printing/send commands for/to R\n\n";
				
				open INSFORR, ">$outdir/$outfileWOpath.$columnheader.$cutofftype.InsFor.R" or die "Can't open '$outfileWOpath.$columnheader.$cutofftype.InsFor.R'\n";
				print INSFORR "##Computing fisher.test for classes ocurring in $columnheader---$cutofftype\n";
				print INSFORR "mat<-read.table(\"$outdir/$outfileWOpath.$columnheader.$cutofftype.matforR.mat\",header=TRUE,row.names=1)
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
	                                       write.table(TableForOutfile.df\$p.val.uncorrected,\"$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalUncorrected\",sep=\"\\t\",quote=FALSE,col.names=NA,row.names=TRUE)
	                                       write.table(p.val.corrected,\"$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalCorrected\",sep=\"\\t\",quote=FALSE,col.names=NA,row.names=TRUE)
					       q()\n";
				close INSFORR;
					       
				###########################
				## Sending to R
				system "R --no-save < $outdir/$outfileWOpath.$columnheader.$cutofftype.InsFor.R";
	
				###########################
				## Obtaining P-values and reformating matrices
				
				print "Obtaining P-values and reformating matrices\n\n";
				
				open UNCORRECTED, "<$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalUncorrected" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalUncorrected'\n";
				$linesHyperUnCorrected = 0;
					while ($line = <UNCORRECTED>) {
						if ($line =~ /^(\S+)(\t)(\S+)/) {
						$Puncorrected = $3;
						$linesHyperUnCorrected++;
						$classid = $hashClassesSentToRPositionName{$linesHyperUnCorrected};
						$hashFinalData{$columnheader}{$cutofftype}{$classid}{puncorrected} = "_$Puncorrected";
						}
					}
				close UNCORRECTED;
					unless ($linesHyperUnCorrected == $classesSentToR) {
					die "ERROR!!! $classesSentToR classes were sent to R for ORA, but $linesHyperUnCorrected uncorrected P-values were retrieved\n";
					}
				
				open CORRECTED, "<$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalCorrected" or die "Can't open '$outdir/$outfileWOpath.$columnheader.$cutofftype.PvalCorrected'\n";
				$linesHyperCorrected = 0;
					while ($line = <CORRECTED>) {
						if ($line =~ /^(\S+)(\t)(\S+)/) {
						$Pcorrected = $3;
						$linesHyperCorrected++;
						$classid = $hashClassesSentToRPositionName{$linesHyperCorrected};
						$hashFinalData{$columnheader}{$cutofftype}{$classid}{pcorrected} = "_$Pcorrected";
						}
					}
				close CORRECTED;
					unless ($linesHyperCorrected == $classesSentToR) {
					die "ERROR!!! $classesSentToR classes were sent to R for ORA, but $linesHyperCorrected corrected P-values were retrieved\n";
					}
					
				system "rm $outdir/$outfileWOpath.$columnheader.$cutofftype.matforR.mat $outdir/$outfileWOpath.$columnheader.$cutofftype.PvalUncorrected $outdir/$outfileWOpath.$columnheader.$cutofftype.PvalCorrected $outdir/$outfileWOpath.$columnheader.$cutofftype.InsFor.R";
				
				}
			}
		}
	}
}

### This will be used to get HYPERG computing time
($year_hyperg2,$month_hyperg2,$day_hyperg2,$hour_hyperg2,$minute_hyperg2,$second_hyperg2) = split ("_", `date +%Y_%m_%d_%H_%M_%S`);


############################
#### Step 5 -- Generating final report
############################
print "Generating final report\n\n";

open OUTMATPUNCORRALL, ">$outdir/$outfileWOpath.ORA.PvalUncorrected.All.mat"                                                     or die "Can't open '$outdir/$outfileWOpath.ORA.PvalUncorrected.All.mat'\n";
open OUTMATPUNCORRFIL, ">$outdir/$outfileWOpath.ORA.PvalUncorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.mat" or die "Can't open '$outdir/$outfileWOpath.ORA.PvalUncorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.mat'\n";
open OUTMATPCORRALL,   ">$outdir/$outfileWOpath.ORA.PvalCorrected.All.mat"                                                       or die "Can't open '$outdir/$outfileWOpath.ORA.PvalCorrected.All.mat'\n";
open OUTMATPCORRFIL,   ">$outdir/$outfileWOpath.ORA.PvalCorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.mat"   or die "Can't open '$outdir/$outfileWOpath.ORA.PvalCorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.mat'\n";

## Summary outfiles
open OUTTABPCORRFIL,       ">$outdir/$outfileWOpath.ORA.PvalCorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.tab"        or die "Can't open '$outdir/$outfileWOpath.ORA.PvalCorrected.$hashParameters{use_ora}$hashParameters{cutoff_ora}.tab'\n";

print OUTMATPUNCORRALL "CutoffType---ClassID---Description";
print OUTMATPUNCORRFIL "CutoffType---ClassID---Description";
print OUTMATPCORRALL   "CutoffType---ClassID---Description";
print OUTMATPCORRFIL   "CutoffType---ClassID---Description";

print OUTTABPCORRFIL       "CutoffType\tClassID\tClassDescription\tColumnName\tPopulation_ClassID_InColumn\tPopulation_ClassID_InCutoff\tP_value\tP_value_corrected\tFold_ratio\n";
print OUTEAPCORRFIL        "ID\tPopulation_ClassID_InColumn\tPopulation_ClassID_InCutoff\tP_value\tP_value_corrected\tFold_ratio\n";
print OUTARRAYSNAPCORRFIL  "ID\tDescription";
print OUTCLASSESNAPCORRFIL "ID\tDescription";

foreach $cutofftype (@CutoffTypes) {
	unless ($hashParameters{$cutofftype} =~ /^na$/i) {
	print OUTARRAYSNAPCORRFIL  "\tPopulation_$cutofftype\_$hashParameters{$cutofftype}";
	print OUTCLASSESNAPCORRFIL "\tPopulation_$cutofftype\_$hashParameters{$cutofftype}";
	}
}
print OUTARRAYSNAPCORRFIL  "\n";
print OUTCLASSESNAPCORRFIL "\n";


#### Printing headers for matrices
foreach $c (1..$NumberOfColumnHeaders) {
$columnheader = $hashHeadersColumns{$c};
	## here restricting to requested columns
	if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {
	
	## here printing column/array headers
	print OUTMATPUNCORRALL "\t$columnheader";
	print OUTMATPUNCORRFIL "\t$columnheader";
	print OUTMATPCORRALL   "\t$columnheader";
	print OUTMATPCORRFIL   "\t$columnheader";
	}
}
print OUTMATPUNCORRALL "\n";
print OUTMATPUNCORRFIL "\n";
print OUTMATPCORRALL   "\n";
print OUTMATPCORRFIL   "\n";

foreach $cutofftype (@CutoffTypes) {
	unless ($hashParameters{$cutofftype} =~ /^na$/i) {
			
		foreach $classid (sort keys %hashClassesSentToRName) {
		
			### Obtain class description
			if ($HashClassRename{$classid}) {
			$classname = $HashClassRename{$classid};
			}else{
			$classname =  "NA";
			}

			### Obtain class maximum population in columns of this matrix
			if ($hashMaxPopulationEachClass{$classid}) {
			$maxpopulation = $hashMaxPopulationEachClass{$classid};
			}else{
			$maxpopulation =  "NA";
			}
		
		print OUTMATPUNCORRALL "$cutofftype---$classid---$classname";
		print OUTMATPUNCORRFIL "$cutofftype---$classid---$classname";
		print OUTMATPCORRALL   "$cutofftype---$classid---$classname";
		print OUTMATPCORRFIL   "$cutofftype---$classid---$classname";
		
			foreach $c (1..$NumberOfColumnHeaders) {
			$columnheader = $hashHeadersColumns{$c};
			## here restricting to requested columns
				if ($hash_numb_cols{$c} or $hashParameters{numb_cols} =~ /^all$/i) {
					if ($hashFinalData{$columnheader}{$cutofftype}{$classid}{puncorrected}) {
					$puncorrected = $hashFinalData{$columnheader}{$cutofftype}{$classid}{puncorrected};
					$puncorrected =~ s/_//;
					}else{
					$puncorrected = "NA";
					}
					if ($hashFinalData{$columnheader}{$cutofftype}{$classid}{pcorrected}) {
					$pcorrected = $hashFinalData{$columnheader}{$cutofftype}{$classid}{pcorrected};
					$pcorrected =~ s/_//;
					}else{
					$pcorrected = "NA";
					}
				
				### Determining if passes the P-value cutoff
				$passesPvalCutoff = 0;
					if ($puncorrected =~ /^na$/i) {
					print OUTMATPUNCORRALL "\t$puncorrected";
					print OUTMATPUNCORRFIL "\t$puncorrected";
					print OUTMATPCORRALL   "\t$pcorrected";
					print OUTMATPCORRFIL   "\t$pcorrected";
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
					
						if ($passesPvalCutoff == 1) {

						$combo = "$columnheader\t$cutofftype\t$classid";
						$hashAllCombosPassingCutoff{$combo} = 1;

						print OUTMATPUNCORRALL "\t$puncorrected";
						print OUTMATPUNCORRFIL "\t$puncorrected";
						print OUTMATPCORRALL   "\t$pcorrected";
						print OUTMATPCORRFIL   "\t$pcorrected";
						
							if ($hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{m}) {
							$m = $hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{m};
							$m =~ s/_//;
							}else{
							$m = "NA";
							}
							
							if ($hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{a}) {
							$a = $hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{a};
							$a =~ s/_//;
							}else{
							$a = "NA";
							}
							
							if ($hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{ratio}) {
							$ratio = $hashDataForHyperG{$columnheader}{$cutofftype}{$classid}{ratio};
							$ratio =~ s/_//;
							}else{
							$ratio = "NA";
							}
						
						print OUTTABPCORRFIL   "$cutofftype\t$classid\t$classname\t$columnheader\t$m\t$a\t$puncorrected\t$pcorrected\t$ratio\n";
						print OUTSIFPCORRFIL   "$classid\t$cutofftype\_$hashParameters{$cutofftype}\t$columnheader\n";
						print OUTEAPCORRFIL    "$classid ($cutofftype\_$hashParameters{$cutofftype}) $columnheader\t$m\t$a\t$puncorrected\t$pcorrected\t$ratio\n";
						
						}else{
						print OUTMATPUNCORRALL "\t$puncorrected";
						print OUTMATPUNCORRFIL "\tNS";
						print OUTMATPCORRALL   "\t$pcorrected";
						print OUTMATPCORRFIL   "\tNS";
						}
					}
				}
			}
		print OUTMATPUNCORRALL "\n";
		print OUTMATPUNCORRFIL "\n";
		print OUTMATPCORRALL   "\n";
		print OUTMATPCORRFIL   "\n";
		}
	}
}
close OUTMATPUNCORRALL;
close OUTMATPUNCORRFIL;
close OUTMATPCORRALL;
close OUTMATPCORRFIL;
close OUTTABPCORRFIL;

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
print OUTRUNTIME "hyperg\t$all_hyperg_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";
close OUTRUNTIME;

print "Took time:
hyperg\t$all_hyperg_seconds\tsecs
overall\t$all_overall_seconds\tsecs\n";



print "\n  Done!!!\n\nCheck outfiles '$outdir/$outfileWOpath.ORA.Pval*'\n
Note in *.na files the population was set with negative for classes and positive for arrays. This is to allow  distinguish as two attributes in Cytoscape\n\n";

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
'infile_matrix' => 1,
'cutoff_neg' => 1,
'cutoff_pos' => 1,
'use_values_or_ranks' => 1,
'infile_classes' => 1,
'restrict_classes' => 1,
'set_min' => 1,
'set_max' => 1,
'use_ora' => 1,
'cutoff_ora' => 1,
'p_correction_test' => 1,
'numb_cols' => 1,
'prefix_outfiles' => 1,
'use_universe' => 1,
);

#######################
#### Starts -- Evaluate parameters

&LoadParameters::Parameters::MainSubParameters(\%hashParametersTolookFor,\@arrayInputtedOneLineCommands);
$Parameters .= "$MoreParameters";

if ($hashParameters{prefix_outfiles} =~ /^default$/) {

## Defining prefix string for OUTFILE
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

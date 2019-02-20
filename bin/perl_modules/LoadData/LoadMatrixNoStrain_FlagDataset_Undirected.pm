############################
## This module populates
## '%hashData{$pair}{$flagDataset} = $score' contains the scores of the matrix for each $pair column/row in $flagDataset
##     In case row or column names repeat in the matrix this module will kill the run
## '%hashHeadersColumnsPerDataset{$c}{$flagDataset}' contains the column headers using numbers from 1..to..N where N is the number of columns in $flagDataset
## '%hashHeadersRowsPerDataset{$c}{$flagDataset}' contains the row headers using numbers from 1..to..M where M is the number of rows in $flagDataset
## '$NumberOfColumnHeaders{$flagDataset}' contains the number of columns in $flagDataset
## '$NumberOfRowHeaders{$flagDataset}'  contains the number of rows in $flagDataset
## '$hashAllPairsIndexed{$pair}' contains ALL pairs of columns/row headers, indexed by a program calling this module (even if the belong to a different '$flaDataset')
##     This may be used as the union of all pairs passed to this module
##
## Note the matrix may contain redundant entries (e.g Gene1[columns] -> Gene2[rows], and Gene2[columns] -> Gene1[rows])
## but only the last one appearing while loading will be indexed. So a symetric matrix is expected
##
## Given the prior note you need to take care when comparing the pairs indexed from more than one database loking at reciprocal pairs as well
##
## The expected format is like:
## ANYTHING   Gene_1    Gene_2    Gene_3
## Gene_1      0.01      0.01      0.03
## Gene_2     -0.01      0.01     -0.15
##
## Note any row/column header may be used. Thus column and row headers may be different e.g. Genes (x-axis) vs. Arrays (y-axis)
##
## The infile may be plain-text either uncompressed or bzipped
##
############################

package LoadData::LoadMatrixNoStrain_FlagDataset_Undirected;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( %hashData %NumberOfColumnHeaders %NumberOfRowHeaders %hashAllPairsIndexed %hashAllHeadersColumns %hashAllHeadersRows );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadMatrixNoStrain_FlagDataset_Undirected {

my($infile_matrix,$flagDataset) = @_;

print "\nLoading values from $flagDataset:\n$infile_matrix\n";

unless (-f $infile_matrix) {
die "Couldn't find infile '$infile_matrix'\n";
}
unless ($flagDataset =~ /\S/) {
die "Couldn't determine dataset flag from '$flagDataset'\n";
}

if ($infile_matrix =~ /\.bz2$/) {
open (INFILE, "bzip2 -qdc $infile_matrix |") or die "Can't open \"$infile_matrix\" (-infile_matrix bzipped)\n";
$dimInfile = `bzcat $infile_matrix | wc`;
}elsif ($infile_matrix =~ /(\.gz$)/) {
open (INFILE, "gzip -qdc $infile_matrix |") or die "Can't open \"$infile_matrix\" (-infile_matrix gzipped)\n";
$dimInfile = `zcat $infile_matrix | wc`;
}else{
open (INFILE, "<$infile_matrix") or die "Can't open \"$infile_matrix\" (-infile_matrix)\n";
$dimInfile = `wc $infile_matrix`;
}

if ($dimInfile =~ /^(\s+)(\d+)/) {
$nrowsInInfile = $2;
	if ($nrowsInInfile > 0) {
	}else{
	die "Couldn't get number of rows for '$infile_matrix'. Does it contain data at the time of this run?\n";
	}
}

$van = 0;
$scoresIndexed = 0;
$scoresIgnored = 0;
$NumberOfColumnHeaders{$flagDataset} = 0;
$NumberOfRowHeaders{$flagDataset} = 0;

while ($line = <INFILE>) {
chomp $line;
$van++;
$c = 0;
	if ($van == 1) {
		if ($line =~ /^(\S+)(\t)(\S.*)/) {
		@ColumnsHeaders = split ("\t", $3);
		}elsif ($line =~ /^(\t)(\S.*)/) {
		@ColumnsHeaders = split ("\t", $2);
		}
		foreach $columnheader (@ColumnsHeaders) {
		$c++;
		$NumberOfColumnHeaders{$flagDataset} += 1;
			if ($hashAlreadyColHeaders{$columnheader}{$flagDataset}) {
			die "Column header '$columnheader' appears more than once in '$flagDataset'. Rename it and try again\n";
			}else{
			$hashAlreadyColHeaders{$columnheader}{$flagDataset} = 1;
			}
		$hashHeadersColumnsPerDataset{$c}{$flagDataset} = $columnheader;
		$hashAllHeadersColumns{$columnheader} = 1;
		}
	}elsif ($line =~ /^(\S+)(\t)(\S.*)/) {
	$rowheader = $1;
	$data = $3;
	@Data = split ("\t", $data);
	$NumberOfRowHeaders{$flagDataset} += 1;
	$NumberOfRowHeadersInThisDataset++;
	$hashHeadersRowsPerDataset{$NumberOfRowHeadersInThisDataset}{$flagDataset} = "$rowheader";
	$hashAllHeadersRows{$rowheader} = 1;

		if ($hashAlreadyRowHeaders{$rowheader}{$flagDataset}) {
		die "Row header '$rowheader' appears more than once in '$flagDataset'. Rename it and try again\n";
		}else{
		$hashAlreadyRowHeaders{$rowheader}{$flagDataset} = 1;
		}
		foreach $score (@Data) {
		$c++;
			
			if ($hashHeadersColumnsPerDataset{$c}{$flagDataset}) {
			$columnheader = $hashHeadersColumnsPerDataset{$c}{$flagDataset};
			
			$pair    = "$rowheader\t$columnheader";
			$pairInv = "$columnheader\t$rowheader";
			
				if ($hashAllPairsIndexed{$pair}) {
				$pairToIndex = "$pair";
				}elsif ($hashAllPairsIndexed{$pairInv}) {
				$pairToIndex = "$pairInv";
				}else{
				$hashAllPairsIndexed{$pair} = 1;
				$pairToIndex = "$pair";
				}
			$hashData{$pairToIndex}{$flagDataset} = "_$score";
			$scoresIndexed++;
			
			}else{
			die "ERROR \"$c\" was not found in hashHeadersColumns\n";
			}
		}
		
		unless ($c == $NumberOfColumnHeaders{$flagDataset}) {
		die "\nERROR!!! the number of column headers ($NumberOfColumnHeaders{$flagDataset}) and data in row $NumberOfRowHeaders{$flagDataset} ($c) is not equal\n\n";
		}

	}else{
	die "ERROR!!! unexpected format in file '$infile_matrix'\n'$line'\n";
	}
}
close INFILE;

print "Done reading matrix ($NumberOfRowHeaders{$flagDataset} rows x $NumberOfColumnHeaders{$flagDataset} columns), $scoresIndexed scores were indexed (including 'NA' values, if any)\n";

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

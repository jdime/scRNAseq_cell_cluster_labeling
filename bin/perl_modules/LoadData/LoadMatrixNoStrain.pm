#########################################
## This module populates
## '%hashData{$columnheader}{$rowheader} = $score' contains the scores of the matrix for each column/row pair
##     with values from a matrix with non redundant neither row nor column names
##     In case row and/or column names are redundat this module will kill the run
## '%hashHeadersColumns{$c}' contains the column headers using numbers from 1..to..N where N is the number of columns
## '%hashHeadersRows{$c}' contains the row headers using numbers from 1..to..M where M is the number of rows
## '$NumberOfColumnHeaders' contains the number of columns
## '$NumberOfRowHeaders'  contains the number of rows
##
## The expected format is like:
## ANYTHING   Array_1   Array_2   Array_3
## Gene_1      0.01      0.01      0.03
## Gene_2     -0.01      0.01     -0.15
##
## Note any row/column header may be used. Thus 'Arrays' (x-axis) may be also genes in a genes vs. genes matrix
##
## The infile may be plain-text either uncompressed or bzipped
##
## Options $capitalizeRows and $capitalizeColumns allow to capitalize row and/or column headers
##
#########################################

package LoadData::LoadMatrixNoStrain;
require Exporter;
require AutoLoader;

use LoadData::HitsToPrintStatus;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( %hashData %hashHeadersColumns %hashHeadersRows $NumberOfColumnHeaders $NumberOfRowHeaders );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadMatrixNoStrain {

my($infile_matrix,$capitalizeRows,$capitalizeColumns) = @_;

%hashData = "";
%hashHeadersColumns = "";
%hashHeadersRows = "";
$NumberOfColumnHeaders = "";
$NumberOfRowHeaders = "";
%hashAlreadyColHeaders = "";
%hashAlreadyRowHeaders = "";

if ($infile_matrix =~ /\.bz2$/) {
open (INFILE, "bzip2 -qdc $infile_matrix |") or die "Can't open \"$infile_matrix\" (-infile_matrix bzipped)\n";
$dimInfile = `bzcat $infile_matrix | wc`;
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

print "\nLoading '$nrowsInInfile' rows from from matrix:\n$infile_matrix\n";

$van = 0;
$scoresIndexed = 0;
$NumberOfColumnHeaders = 0;
$NumberOfRowHeaders = 0;
while ($line = <INFILE>) {
chomp $line;
$van++;
$c = 0;

 	if ($HitsToPrintStatus{$van}) {
 	print "$van lines readed\n";
 	}
	
	if ($van == 1) {
		if ($capitalizeColumns =~ /^y$/i) {
		$line =~ tr/[a-z]/[A-Z]/;
		}
	$columnHeadersWithoutSpaces = $line;
	$columnHeadersWithoutSpaces =~ s/ +/_/;

		if ($columnHeadersWithoutSpaces =~ /^(\S* *\S+)(\t)(\S.*)/) {
		@ColumnsHeaders = split ("\t", $3);
		}elsif ($columnHeadersWithoutSpaces =~ /^(\t)(\S.*)/) {
		@ColumnsHeaders = split ("\t", $2);
		}else{
		die "\nERROR!!! couldn't identify column headers\n";
		}
		
		foreach $columnheader (@ColumnsHeaders) {
		$c++;
		$NumberOfColumnHeaders++;
		$columnheader =~ s/\(\d+\)//g;

			if ($hashAlreadyColHeaders{$columnheader}) {
			die "Column header '$columnheader' appears more than once. Rename it and try again\n";
			}else{
			$hashAlreadyColHeaders{$columnheader} = 1;
			}
		$hashHeadersColumns{$c} = $columnheader;
		}
	}elsif ($line =~ /^(\S+)(\t)(\S.*)/) {
	$rowheader = $1;
	$data = $3;
	
		if ($capitalizeRows =~ /^y$/i) {
		$rowheader =~ tr/[a-z]/[A-Z]/;
		}

	@Data = split ("\t", $data);
	$NumberOfRowHeaders++;
	$hashHeadersRows{$NumberOfRowHeaders} = "$rowheader";

		if ($hashAlreadyRowHeaders{$rowheader}) {
		die "Row header '$rowheader' appears more than once. Rename it and try again\n";
		}else{
		$hashAlreadyRowHeaders{$rowheader} = 1;
		}
		foreach $score (@Data) {
		$c++;
			if ($hashHeadersColumns{$c}) {
			$columnheader = $hashHeadersColumns{$c};
			$hashData{$columnheader}{$rowheader} = "_$score";
			$scoresIndexed++;
			}else{
			die "ERROR \"$c\" was not found in hashHeadersColumns\n";
			}
		}

		unless ($c == $NumberOfColumnHeaders) {
		die "\nERROR!!! the number of column headers ($NumberOfColumnHeaders) and data in row $NumberOfRowHeaders ($c) is not equal\n\n";
		}

	}else{
	die "ERROR!!! unexpected format in file '$infile_matrix'\n'$line'\n";
	}
}
close INFILE;

print "Done reading matrix ($NumberOfRowHeaders rows x $NumberOfColumnHeaders columns), $scoresIndexed scores were indexed (including 'NA' values, if any)\n";

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

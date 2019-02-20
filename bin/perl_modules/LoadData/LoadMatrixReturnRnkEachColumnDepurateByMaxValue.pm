############################
#### sub(LoadMatrixReturnRnkEachColumnDepurateByMaxValue) indexes scores for pairs of identifiers, e.g. arrays (x-axis) and genes (y-axis); or genes vs. genes (in both axes)
############################
###
### Expected format:
### COMPLETE   Array_1   Array_2   Array_3
### Gene_1     0.01      0.01      0.03
### Gene_2:A   -0.01     0.01     -0.15
### Gene_2:B   -0.03     0.04      0.05
###
### The string COMPLETE may exist or not, or may be any other string
###
### NOTES:
###       1) Row headers will be capitalized
###       2) Row headers with non-alpha numeric characters will be replaced by low dash '_'. If that happens $ChangedColNames will be > 0, and $outfileOriginalColNames will store the name of the file with old=>new headers
###          This helps for example to avoid crashes with R
###       3) The infile_matrix may be plain-text either uncompressed or gunzipped or bzipped
###       4) Redundant column headers are not allowed
###
### Redundant pairs ColumnName/RowName will be depurated retaining only the higher ab(value). For example:
### Gene_2 appears twice. The ':' are used as separators and everything after the ':' will be removed from RowName. Thus the final scores for Gene_2 pairs will be:
### Array_1 vs. Gene2 => -0.03
### Array_2 vs. Gene2 =>  0.04
### Array_3 vs. Gene2 => -0.15
###
############################
###
### Returns a list of Genes and scores for each Array, in a Array*.rnk format, e.g. useful for GSEA, e.g.:
### Gene_1   0.01
### Gene_2  -0.01
###
############################
### The infile may be a gzipped, bzipped or plain text file
############################

package LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue;
require Exporter;
require AutoLoader;

use ReformatNumbers::NAvalues;
use ReformatPerlEntities::ChangeCharactersByLowDash;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( LoadMatrixReturnRnkEachColumn %hashAllRowNames %hashAllColumnNames %hashRowsNames %hashColumnsNames $NumberOfRowHeaders $NumberOfColumnHeaders %hashAllReplicatesLoaded %hashCountScoresPassingCutoffsWONas $ChangedColNames $outfileOriginalColNames );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadMatrixReturnRnkEachColumnDepurateByMaxValue {

($path_outfiles,$infile_matrix,$cutoff_from_input,$abs_from_input,$more_or_less,$use_universe,$row_names_for_intersection) = @_;

ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($infile_matrix);

@CheckDefined = ($path_outfiles,$infile_matrix,$cutoff_from_input,$abs_from_input,$more_or_less,$use_universe,$row_names_for_intersection);
$c = 0;
foreach $to_check (@CheckDefined) {
$c++;
	if ($to_check =~ /\S/) {
	### Do nothing
	}else{
	die "\nERROR!!! in module LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue the '$c' value inputted by \@_ is empty\n";
	}
}

`mkdir -p $path_outfiles`;

(my %hashData = "");

if ($use_universe =~ /^i$/i) {
%hash_row_names_for_intersection = %{$row_names_for_intersection};
}elsif  ($use_universe =~ /^n1$/i) {
### Do nothing
}else{
die "\n\nERROR!!! unexpected option \$use_universe '$use_universe'\n\n";
}


### Hashing matrix
print "Reading Matrix from $infile_matrix with LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue\n";

$van = 0;
$ChangedColNames = 0;
$NumberOfRowHeaders = 0;
$NumberOfColumnHeaders = 0;

if ($infile_matrix =~ /\.bz2$/) {
open (MATRIX, "bzip2 -qdc $infile_matrix |") or die "Can't open \"$infile_matrix\" (-infile_matrix bzipped)\n";
}elsif ($infile_matrix =~ /\.gz$/ or $infile_matrix =~ /\.Z$/) {
open (MATRIX, "gzip -qdc $infile_matrix |") or die "Can't open \"$infile_matrix\" (-infile_matrix gunzipped)\n";
}else{
open (MATRIX, "<$infile_matrix") or die "Can't open \"$infile_matrix\" (-infile_matrix)\n";
}

while ($line = <MATRIX>) {
chomp $line;
$van++;
$c = 0;
$r++;


	if ($line =~ /^(\S+)(\t)(\S.*)/) {
	@Data = split ("\t", $3);
	$rowname = $1;
	$rowname =~ tr/[a-z]/[A-Z]/;
	}elsif ($line =~ /^(\t)(\S.*)/) {
	@Data = split ("\t", $2);
	$rowname = "NA";
	}

	if ($van == 1) {
		foreach $colname (@Data) {
		$c++;
		
		$NumberOfColumnHeaders++;
		$OriginalColName = $colname;
		$colname = &ReformatPerlEntities::ChangeCharactersByLowDash::ChangeCharactersByLowDash($colname);
			if ($hashAllColumnNames{$colname}) {
			die "\nERROR!!! column header '$colname' appears more than once. Duplicted headers are not allowed\n";
			}else{
			$ColumnName = "$colname";
			}
			
		$hashColumnsNames{$NumberOfColumnHeaders} = "$ColumnName";
		$hashAllColumnNames{$ColumnName} = 1;
		$hashOriginalToNewColumnNames{$NumberOfColumnHeaders} = "$OriginalColName\t$ColumnName";
		
			unless ($colname eq $OriginalColName) {
			$ChangedColNames++;
			}
		}
		
	### Generating an outfile if ColumnNames were modfied
		if ($ChangedColNames > 0 ) {
		$outfileOriginalColNames = "$path_outfiles/$outfileWOpath.OriginalColNames";
		open ORIGINALCOLNAMES, ">$outfileOriginalColNames" or die "Can't open '$outfileOriginalColNames'\n";
		print ORIGINALCOLNAMES "# '$ChangedColNames' column names changed\n#Original_ColumnName\tNew_ColumnName\n";
			foreach $colnumber (1..$NumberOfColumnHeaders) {
			print ORIGINALCOLNAMES "$hashOriginalToNewColumnNames{$colnumber}\n";
			}
		close ORIGINALCOLNAMES;
		}
	
	}else{

		if (($use_universe =~ /^i$/i && $hash_row_names_for_intersection{$rowname}) or ($use_universe =~ /^n1$/i)) {
		
		
			### Here just printing warnings about duplicate rownames. But the actuall indexing is below and it will depend on the values of identical pairs columname/rowname
			if ($hashAllRowNames{$rowname}) {
			print "WARNING!!! row name '$rowname' appears more than once. Only the max. abs(value) for each pair ColumnName/RowName will be retained\n";
			}else{
			$hashAllRowNames{$rowname} = 1;
			}
		
			$NumberOfRowHeaders++;
			$hashRowsNames{$NumberOfRowHeaders} = $rowname;
		
			foreach $score (@Data) {
			$score =~ tr/[a-z]/[A-Z]/;
			$c++;
			$ColumnName = $hashColumnsNames{$c};
			$passesCutoffFromInput = 0;
			
				### Filtering out 'NA' and empty values
				unless ($hashNAvalues{$score} or $score =~ /^$/) {
					
					### Filtering out by input score (if applicable)
					if ($cutoff_from_input =~ /(^na$)|(^none$)/i) {
					$passesCutoffFromInput = 1;
					}else{
						if ($abs_from_input =~ /^y$/i) {
							if ($more_or_less =~ /more/i) {
								if (abs($score) >= $cutoff_from_input) {
								$passesCutoffFromInput = 1;
								}
							}elsif ($more_or_less =~ /less/i) {
								if (abs($score) <= $cutoff_from_input) {
								$passesCutoffFromInput = 1;
								}
							}else{
							die "\nERROR!!! in module LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue comparison couldn't be made, either '>=' or '<=', defined by 'more' or 'less' was expected, but '$more_or_less' was passed instead\n\n";
							}
						}else{
							if ($more_or_less =~ /more/i) {
								if ($score >= $cutoff_from_input) {
								$passesCutoffFromInput = 1;
								}
							}elsif ($more_or_less =~ /less/i) {
								if ($score <= $cutoff_from_input) {
								$passesCutoffFromInput = 1;
								}
							}else{
							die "\nERROR!!! in module LoadData::LoadMatrixReturnRnkEachColumnDepurateByMaxValue comparison couldn't be made, either '>=' or '<=', defined by 'more' or 'less' was expected, but '$more_or_less' was passed instead\n\n";
							}
						}
					}
					
					if ($passesCutoffFromInput == 1) {
						if ($hashData{$ColumnName}{$rowname}) {
						$oldValue = $hashData{$ColumnName}{$rowname};
						$oldValue =~ s/^_//;
							if (abs($oldValue) < abs($score)) {
							$hashData{$ColumnName}{$rowname} = "_$score";
							}
						}else{
						$hashData{$ColumnName}{$rowname} = "_$score";
						}
					}
				}
			}
			
			unless ($c == $NumberOfColumnHeaders) {
			die "\nERROR!!! the number of column headers ($NumberOfColumnHeaders) and data in row $NumberOfRowHeaders ($c) is not equal\n\n";
			}
		}
	}
}
close MATRIX;

#### Generating outfiles for each column
foreach $ColumnName (keys %hashData) {
	if ($ColumnName =~ /\S/) {
	print "Write '$path_outfiles/$ColumnName.rnk\n";
	open RANK_OUT, ">$path_outfiles/$ColumnName.rnk" or die "Can't open '$path_outfiles/$ColumnName.rnk'\n";
		foreach $rowname (keys %{$hashData{$ColumnName}}) {
		$score = $hashData{$ColumnName}{$rowname};
		$score =~ s/^_//;
		print RANK_OUT "$rowname\t$score\n";
		}
	close RANK_OUT;
	}
}


}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

#########################
### Subroutine 'ObtainOutfileWOpath'
### takes as input a string and removes both the prefix (e.g. path to a file)
### and a variety of suffix strings [e.g. file extension(s)] specified at %hashSuffixToRemove
### to return a $outfileWOpath which may be used as a template for outfile names
#########################

package ReformatPerlEntities::ObtainOutfileWOpath;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( ObtainOutfileWOpath $outfileWOpath $pathToOutfile );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub ObtainOutfileWOpath {

### Put here all extensions that want to be removed from inputted string
%hashSuffixToRemove = (
'na' => 1,
'ea' => 1,
'all' => 1,
'fnu' => 1,
'fna' => 1,
'faa' => 1,
'fas' => 1,
'fsa' => 1,
'fasta' => 1,
'fastq' => 1,
'list' => 1,
'mat' => 1,
'net' => 1,
'obo' => 1,
'tab' => 1,
'tmp' => 1,
'txt' => 1,
'xlsx' => 1,
'clstr' => 1,
'ptt' => 1,
'reformated' => 1,
'jpg' => 1,
'jpeg' => 1,
'png' => 1,
'bmp' => 1,
'csv' => 1,
'tsv' => 1,
'mtx' => 1,
'gmt' => 1,
'gtf' => 1,
'gff3' => 1,
'gff' => 1,
);

my($inString) = @_;
$outfile = $inString;
	if ($outfile =~ /^(\S+\/)/) {
	$pathToOutfile = $1;
	}else{
	$pathToOutfile = "";
	}
$pathToOutfile =~ s/\/$//;

$outfile =~ s/\S+\///; ## removing any path/directories
$outfile =~ s/(\.bz2|\.gz|\.Z)$//;
$outfile =~ s/\.tar$//;


	if ($outfile =~ /^(\S+)(\.)([a-z]+)$/i) {
	$outfilePart1 = $1;
	$outfilePart2 = $3;
	$outfilePart2 =~ tr/[A-Z]/[a-z]/;
		if ($hashSuffixToRemove{$outfilePart2}) {
		$outfileWOpath = $outfilePart1;
		}else{
		$outfileWOpath = "$outfilePart1.$outfilePart2";
		}
	}else{
	$outfileWOpath = $outfile;
	}
return $outfileWOpath;

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

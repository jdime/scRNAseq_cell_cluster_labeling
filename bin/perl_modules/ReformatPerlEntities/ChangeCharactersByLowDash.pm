#########################
### Subroutine 'ChangeCharactersByRequestedCharacter' changes all characters defined
### in $outputString =~ s/.../_/g;
### by another low dash '_'
### Then returns the modified $outputString
#########################

package ReformatPerlEntities::ChangeCharactersByLowDash;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( ChangeCharactersByLowDash ChangeCharactersByLowDash_ForFread $outputString );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub ChangeCharactersByLowDash {

my ($inputString) = @_;
my ($outputString) = "";

$outputString = $inputString;
$outputString =~ s/\s+|;|:|,|-|_|"|'|\(|\)|\[|\]|\{|\}|\/|\\|\||\+|!|&/_/g;
$outputString =~ s/_+/_/g;
	unless ($outputString =~ /^_$/) {
	$outputString =~ s/^_//;
	$outputString =~ s/_$//;
	}

	if ($outputString =~ /\S/) {
	}else{
	die "\nERROR!!! module ReformatPerlEntities::ChangeCharactersByLowDash in string\n'$inputString' replacing characters by '_'\nreturns an empty string\n'$outputString'\n\n";
	}
return "$outputString";
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub ChangeCharactersByLowDash_ForFread {

### Some characters in column headers are not properly handled and couldn't be escaped by R's command 'fread' library("data.table")
### Thus they can be converted to '_'. Dots are OK though.

my ($inputString) = @_;
my ($outputString) = "";

$outputString = $inputString;
$outputString =~ s/ |;|:|,|-|_|"|'|\(|\)|\[|\]|\{|\}|\/|\\|\||\+|!/_/g;
	unless ($outputString =~ /^_$/) {
	$outputString =~ s/^_//;
	$outputString =~ s/_$//;
	}

	if ($outputString =~ /\S/) {
	}else{
	die "\nERROR!!! module ReformatPerlEntities::ChangeCharactersByLowDash_ForFread in string\n'$inputString' replacing characters by '_'\nreturns an empty string\n'$outputString'\n\n";
	}
return $outputString;
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



1;

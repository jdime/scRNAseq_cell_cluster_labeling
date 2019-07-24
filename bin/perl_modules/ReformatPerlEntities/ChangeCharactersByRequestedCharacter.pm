#########################
### Subroutine 'ChangeCharactersByRequestedCharacter' changes all characters defined
### in $outputString =~ s/.../$characterReplacing/g
### by $characterReplacing
### Then returns the modified $outputString
###
### Variable $compressYN allows to replace consecutive $characterReplacing characters
### by a single one
#########################

package ReformatPerlEntities::ChangeCharactersByRequestedCharacter;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( ChangeCharactersByRequestedCharacter $outputString );

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub ChangeCharactersByRequestedCharacter {

my ($inputString,$characterReplacing,$compressYN) = @_;
my ($outputString) = "";

$outputString = $inputString;
$outputString =~ s/\s+|;|:|,|-|_|"|'|\(|\)|\[|\]|\{|\}|\/|\\|\||\+|!|&/$characterReplacing/g;

	if ($compressYN =~ /^y$/i) {
	$outputString =~ s/$characterReplacing+/$characterReplacing/g;
	}elsif ($compressYN =~ /^n$/i) {
	### Do nothing
	}else{
	die "\n\nERROR!!! unexpected compress option '$compressYN' in ReformatPerlEntities::ChangeCharactersByRequestedCharacter\n\n";
	}

	unless ($outputString =~ /^$characterReplacing$/) {
	$outputString =~ s/^$characterReplacing//;
	$outputString =~ s/$characterReplacing$//;
	}

	if ($outputString =~ /\S/) {
	}else{
	die "\nERROR!!! module ReformatPerlEntities::ChangeCharactersByRequestedCharacter in string\n'$inputString' replacing characters by '$characterReplacing'\nreturns an empty string\n'$outputString'\n\n";
	}
return "$outputString";
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

#########################
## This module populates %hashNAvalues which contains strings
## that can be used as 'NA' values (e.g. empty fields)
##
## Note: whenever you call this module make sure you use upper-case
##       strings, e.g. using 'tr' command:
##       $string =~ tr/[a-z]/[A-Z];
#########################

package ReformatNumbers::NAvalues;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw(
IndexesHashNAvalues
%hashNAvalues
$nas
);

&IndexesHashNAvalues;

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub IndexesHashNAvalues {

%DefineNaValuesHere = (
'N/A' => 1,
'NA' => 1,
'NALINKAGE' => 1,
'NASELF' => 1,
'NAMISSING' => 1,
'NADEFECTIVE' => 1,
'NANA' => 1,
'NAEMPTY' => 1,
'NANONREP' => 1,
'NAN' => 1,
'INF' => 1,
'-INF' => 1,
'NODATA' => 1,
'ND' => 1,
'NULL' => 1,
'#DIV/0!' => 1,
'#NAME?' => 1,
'SAMEBARCODE' => 1,
'SAMESTRAIN' => 1,
'SAMEPOSITION' => 1,
'NOTTESTED' => 1,
);

foreach $na (sort keys %DefineNaValuesHere) {
$na =~ tr/[a-z]/[A-Z]/;
$hashNAvalues{$na} = 1;
$nas .= ", \"$na\"";
}

print "\nThe following strings will be considered 'NA' values. Make sure your program uses uppercase strings to match them\n";
$nas =~ s/^, //;
print "$nas\n\n";

return $nas;

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

#########################
#### sub(LoadClasses) indexes the number and identifiers of members of classes in *.gmt format
#### Also provides the name of each class
#### Also generates an outifle with classes matching requested population in parameters '$set_min' and '$set_max'
####      and NOT contained in file $restrict_classes.
####      These two filters may help to save running time
#### 
#### An outfile $infile_classes.Filtered.tmp will be created with classes/genes passing above filters, and also with class identifier
####      and names capitalized, and non alpha numeric symbols. This is neccesary for example to run GSEA and match class descriptions
####
############################
#### Files MUST have ONE CLASS PER LINE, and can include or not class_ID and class_NAME
#### indicated by 'yes' or 'no' by $includes_class_id_and_name_YN
###### e.g. 
## GO:0000103	6.sulfate.assimilation	cysC__b2750__JW2720	cysD__b2752__JW2722
## GO:0001101	5.response.to.acid	adiA__b4117__JW5731	adiC__b4115__JW4076
#########################

package LoadData::LoadClasses;
require Exporter;
require AutoLoader;

use ReformatPerlEntities::ChangeCharactersByRequestedCharacter;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw(
LoadRestrictClasses
LoadClasses
%HashEachKeyClassesPassingCutoff
%hashOrderOfClassesRetained
$ClassesRetained
%HashClassRename
%HashClassKeys
%HashOfHashClassEachKeys
%HashEachKeyClasses
%HashClassData
%HashClassesAllLoaded
%hashClassesRetained
$Extras
$ReformattedClasses
$count_genes_tested_from_gene_list
$count_genes_tested_from_gene_list_and_classes
);

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadClasses {

my($infile_classes,$restrict_classes,$set_min,$set_max,$infile_tested_genes,$column_tested_genes,$includes_class_id_and_name_YN,$outfile_filtered_classes) = @_;

unless ($infile_classes =~ /\S/ && $restrict_classes =~ /\S/ && $set_min =~ /^\d+$/ && $set_max =~ /^\d+$/ && $infile_tested_genes =~ /\S/ && $column_tested_genes =~ /\S/ && $includes_class_id_and_name_YN =~ /^(yes|no)$/i && $outfile_filtered_classes =~ /\S/) {
die "Unexpected format in options for LoadData::LoadClasses:
infile_classes           = '$infile_classes'
restrict_classes         = '$restrict_classes'
set_min                  = '$set_min'
set_max                  = '$set_max'
infile_tested_genes      = '$infile_tested_genes'
column_tested_genes      = '$column_tested_genes'
include_class_name       = '$includes_class_id_and_name_YN'
outfile_filtered_classes = '$outfile_filtered_classes'\n";
}

&LoadRestrictGenes($infile_tested_genes,$column_tested_genes);
&LoadRestrictClasses($restrict_classes);

print "Loading Classes from \"$infile_classes\"\n";
open CLASSES,    "<$infile_classes"           or die "Can't open '$infile_classes' (-infile_classes)\n";
open CLASSESOUT, ">$outfile_filtered_classes" or die "Can't open '$outfile_filtered_classes' (outfile_classes filtered by population) \n";

$ReformattedClasses = "$outfile_filtered_classes";

$classnumber = 0;
$ClassesLoaded = 0;
$ClassesIgnored = 0;
$ClassesFilteredByPopulation = 0;
$ClassesRetained = 0;
$van = 0;
$count_genes_tested_from_gene_list_and_classes = 0;
while ($line = <CLASSES>) {
chomp $line;

	unless ($line =~ /^#/) {
	$line =~ s/\s*\t\s*/\t/g;
	$line =~ tr/[a-z]/[A-Z]/;
	$van++;
	
	$ToEnterIndexing = 0;
	
		if ($includes_class_id_and_name_YN =~ /^yes$/i) {
			if ($line =~ /^(\S+)(\t)(\S+)(\t)(\S.*)/) {
			$classid = $1;
			$classname = $3;
			$classnodes = $5;
			$ToEnterIndexing = 1;
			}
		}elsif ($includes_class_id_and_name_YN =~ /^no$/i) {
			if ($line =~ /^(\S.*)/) {
			$classid = $van;
			$classname = $van;
			$classnodes = $1;
			$ToEnterIndexing = 1;
			}
		}
		
		if ($ToEnterIndexing == 1) {

		$PassesRestriction = 0;

			if ($restrict_classes =~ /^all$/i) {
			$PassesRestriction = 1;
			}else{
				if ($HashClassData{$classid}{restricted}) {
				$PassesRestriction = 1;
				}
			}
			
			if ($PassesRestriction == 1) {
			$ClassesLoaded++;

			$classname =~ s/(<\/*SUP>|<\/*I>)/_/g;
			$classname = ReformatPerlEntities::ChangeCharactersByRequestedCharacter::ChangeCharactersByRequestedCharacter($classname,"_","y");

			$classid =~ s/(<\/*SUP>|<\/*I>)/_/g;
			$classid = ReformatPerlEntities::ChangeCharactersByRequestedCharacter::ChangeCharactersByRequestedCharacter($classid,"_","y");
			
			$classname =~ tr/[a-z]/[A-Z]/;
			$classid =~ tr/[a-z]/[A-Z]/;
			
			$HashClassRename{$classid} = $classname;
			$HashClassKeys{$classid} = "$classnodes";
			@ClassNodes = split ("\t", $classnodes);
			%hashAllNodesThisClass = "";
			
				foreach $key1 (@ClassNodes) {
					if ($key1 =~ /\S/) {
						if ($HashGenesTested{$key1} or $infile_tested_genes =~ /^NA$/i) {
						$hashAllNodesThisClass{$key1} = 1;
	
							unless ($HashOfHashClassEachKeys{$classid}{$key1}) {
							$HashOfHashClassEachKeys{$classid}{$key1} = 1;
							$HashEachKeyClasses{$key1} .= "$classid,";
							$HashClassesAllLoaded{$classid} += 1;
									
								unless ($HashGenesTested{$geneid}{infile_tested_genes_and_classes}) {
								$HashGenesTested{$geneid}{infile_tested_genes_and_classes} = 1;
								$count_genes_tested_from_gene_list_and_classes++;
								}
							}
						}
					}
				}
				
				### Sorting alphanumerically nodes in this class and getting number of nodes per class
				$classNodesSortedAlphanumerically = "";
				foreach $key1 (sort keys %hashAllNodesThisClass) {
					if ($key1 =~ /\S/) {
					$HashClassData{$classid}{nodesAllCount} += 1;
					$HashClassData{$classid}{nodesAllList} .= "$key1\n";
					$classNodesSortedAlphanumerically      .= "\t$key1";
					}
				}
				$classNodesSortedAlphanumerically =~ s/^\t//;

				if (($HashClassData{$classid}{nodesAllCount} >= $set_min) && ($HashClassData{$classid}{nodesAllCount} <= $set_max)) {
				$hashClassesRetained{$classid} = 1;
					if ($includes_class_id_and_name_YN =~ /^yes$/) {
					print CLASSESOUT "$classid\t$classname\t$classNodesSortedAlphanumerically\n";
					@ClassNodesSortedAlphanumerically = split ("\t", $classNodesSortedAlphanumerically);
					
						foreach $key1 (@ClassNodesSortedAlphanumerically) {
						$HashEachKeyClassesPassingCutoff{$key1} .= "$classid,";
						}
					
					}else{
					print CLASSESOUT "$classNodesSortedAlphanumerically\n";
					}
					
				$ClassesRetained++;
				$hashOrderOfClassesRetained{$ClassesRetained} = $classid;
				}else{
				$ClassesFilteredByPopulation++;
				}
			}else{
			$ClassesIgnored++;
			}
		}
	}
}
close CLASSES;
close CLASSESOUT;

print "\n$ClassesLoaded classes were loaded
$ClassesIgnored classes were ignored
$ClassesFilteredByPopulation classes were filtered by population (min=$set_min, max=$set_max)
$ClassesRetained classes were retained and written to outfile '$infile_classes.Filtered.gmt'\n\n";

$Extras .= "$ClassesLoaded classes were loaded
$ClassesIgnored classes were ignored
$ClassesFilteredByPopulation classes were filtered by population (min=$set_min, max=$set_max)
$ClassesRetained classes were retained and written to outfile '$infile_classes.Filtered.gmt'\n";
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadRestrictClasses {
############################
#### Loading restricted classes
############################
## This may be used to filter the list of classes inputted from -infile_classes. E.g. when using GO terms comparing several pairs of big classes (even using cutoffs)
## causes the program crashes by exhausting CPU memory
###### The infile is a list of class identifier that must match with -infile_classes (first column), e.g.
## GO:0000103
## GO:0001101
############################

my($restrict_classes) = @_;

if ($restrict_classes =~ /^all$/i) {
print "All classes from -infile_classes will be passed to next filtering (by population size)\n";
}else{

	if ($includes_class_id_and_name_YN =~ /^no$/i) {
	die "\n\nERROR!!! uncompatible options restrict_classes = '$restrict_classes' and includes_class_id_and_name_YN = '$includes_class_id_and_name_YN'\n";
	}

	print "Loading Restricting classes from \"$restrict_classes\"\n";
	open RESTRICTCLASSES, "<$restrict_classes" or die "Can't open \"$restrict_classes\" (infile -restrict_classes)\n";
	$restrictedclassnumber = 0;
	while ($line = <RESTRICTCLASSES>) {
	chomp $line;
	$line =~ tr/[a-z]/[A-Z]/;
		unless ($line =~ /^#/) {
			if ($line =~ /^(\S+)/) {
			$classid = $1;
			$HashClassData{$classid}{restricted} = 1;
			$restrictedclassnumber++;
			}
		}
	}
	close RESTRICTCLASSES;
	print "$restrictedclassnumber classes were loaded. Those will be indexed from -infile_classes and will be passed to next filtering (by population)\n";
}

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub LoadRestrictGenes {
############################
#### Loading restricted genes
############################
## This may be used to filter the list of genes inputted from $column_tested_genes of $infile_tested_genes
## when the enrichment statistics should be restricted to those genes
## The infile is a list of gene identifiers, one-per row
############################

my($infile_tested_genes,$column_tested_genes) = @_;

if ($infile_tested_genes =~ /^NA$/i) {
print "All genes from -infile_classes will be passed to next filtering (by population size)\n";
}else{

	print "Loading Restricting genes from '$infile_tested_genes' column '$column_tested_genes'\n";
	open RESTRICTGENES, "<$infile_tested_genes" or die "Can't open '$infile_tested_genes' (infile -infile_tested_genes)\n";
	$count_genes_tested_from_gene_list = 0;
	$van = 0;
	while ($line = <RESTRICTGENES>) {
	chomp $line;
	$line =~ tr/[a-z]/[A-Z]/;
	$van++;
	
		unless ($line =~ /^#/ or ($van == 1 && $line =~ /^ID|^Protein|^Gene|^Node/i)) {
		@arr = split ("\t", $line);
		$column_tested_genes_m1 = $column_tested_genes - 1;
		$geneid = @arr[$column_tested_genes_m1];
			unless ($HashGenesTested{$geneid}{infile_tested_genes}) {
			$HashGenesTested{$geneid}{infile_tested_genes} = 1;
			$count_genes_tested_from_gene_list++;
			}
		}
	}
	close RESTRICTGENES;
	print "$count_genes_tested_from_gene_list genes were loaded. Those will be indexed from -infile_tested_genes and will be passed to next filtering (by population)\n";
}

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1;

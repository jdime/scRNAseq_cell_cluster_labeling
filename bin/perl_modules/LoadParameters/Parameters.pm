
### Modules LoadParameters::Evaluate_Definitions and LoadParameters::Parameters
### work together to request and evaluate parameters for Perl programs

### They were inspired by libraries "argparse" from Python and "optparse" from R

### Module LoadParameters::Parameters sorts the parameter and sends it to the corresponding evaluates that one-command-line parameters provided by the User
### meet the expected format (numeric, integer, strings, <comma> separated options, etc)

package LoadParameters::Parameters;
require Exporter;
require AutoLoader;

use LoadParameters::Evaluate_Definitions;
use PathsDefinition::PathsToInputs;

### NOTE: hashes generated in LoadParameters::Evaluate_Definitions need to be #EXPORTed here too to be used by Perl programs

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw(
MainSubParameters
%hashParameters
$MoreParameters
%hashColumnKeyToColumnNumber
%hash_parameter_values_to_include
%hash_numb_cols
%hash_numb_rows
%hash_outfile_options
%hash_step_numbers_to_run
%hash_numb_rounds
%hash_rank_classes_to_include
%hash_cutoffs_numeric
%hash_cutoffs_negative
%hash_cutoffs_positive
%hash_cutoffs_negative_or_positive
%hashInflations
%hash_distance_measures
%hash_clustering_methods
$TypeOfTipForRobot_Label
$Horizontal_Offset
);

sub MainSubParameters {

#### Loading inputted parameters

my($In_hashParametersTolookFor,$In_arrayInputtedOneLineCommands) = @_;
%hashParametersTolookFor = %{$In_hashParametersTolookFor};   ## here assigning the incoming hash (element 0 from @_);
@InputtedOneLineCommands = @{$In_arrayInputtedOneLineCommands}; ## here assigning the incoming array (element 1 from @_);

#### Indexing parameters
	
	$c = 0;
	$CountMissingParamaters = 0;
	foreach $parameter (@InputtedOneLineCommands) {
	$c++;
	
		if ($parameter =~ /^-/) {
			if ($parameter =~ /^-\d+/) { ## this is because a parameter values like '-aggrv_cutoff -2' might be considered as two parameters because the '-'
			$MoreParameters .= "$parameter\n";
			}else{
			$MoreParameters .= "$parameter \t";
			}
		}else{
		$MoreParameters .= "$parameter\n";
		}
	
		if ($parameter =~ /^(-)(\S+)$/) {
		$parameterKey = $2;
		$parameterValue = @InputtedOneLineCommands[$c];
			if ($parameterValue =~ /(^-*0+$)/) {
			$parameterValue = "0.0";
			}
		$hashParameters{$parameterKey} = "$parameterValue";
		}
	}
	
#### Checking parameters for further details
	
	foreach $parameterKey (sort keys %hashParametersTolookFor) {
		if ($hashParameters{$parameterKey}) {
		$parameterValue = $hashParameters{$parameterKey};
	
		#### These are parameters which need further processing
		#### Also you need to define such details in module LoadParameters::EvaluateDefinitions
		#### Also you need to return any useful hash, array, etc from such other modules by @EXPORT
		
		#### IMPORTANT: keep this row to START parsing of 'make_list_of_perl_dependencies_in_perl_programs.pl' script

			if ($parameterKey =~ /(^infile)/) {
			### these are parameters which only specified path/files do exist
				### This is to allow Evaluate_Definitions::Evaluate_file_exist use '-f' when the input starts with '~/' instead the root/home/Users path
				if ($parameterValue =~ /^~\//) {
				$PathPrefix = "/$Users_home/$DefaultUserName/";
				$parameterValue =~ s/^~\//$PathPrefix/;
				$hashParameters{$parameterKey} = $parameterValue;
				}
			LoadParameters::Evaluate_Definitions::Evaluate_file_exist($parameterKey,$parameterValue);

			}elsif ($parameterKey =~ /(^prefix)/) {
			### these are parameters which only need to check that are defined 
			LoadParameters::Evaluate_Definitions::Evaluate_were_defined($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cex_labels$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^classes_type$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cutoff_neg$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cutoff_ora$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cutoff_pos$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cutoff_print_p$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^cutoff_print_q$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^lwd$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^nperm$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^p_correction_test$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^path_outfiles$)/) {
			### these are parameters which only need specified path/directories do exist
			$hashParameters{$parameterKey} =~ s/\/$//;
			LoadParameters::Evaluate_Definitions::Evaluate_directory_exist($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^print_plot_ticks_labels_legend$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^use_graphics_device$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^use_ora$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^use_universe$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);
			}elsif ($parameterKey =~ /(^use_values_or_ranks$)/) {
			LoadParameters::Evaluate_Definitions::Evaluate_value($parameterKey,$parameterValue);

		#### IMPORTANT: keep this row to END parsing of 'make_list_of_perl_dependencies_in_perl_programs.pl' script

			}else{
			die "\nERROR!!! parameter '$parameterKey' is not definied in module LoadParameters::Parameters\n";
			}

		}else{
		$hashMissingParameters{$parameterKey} = 1;
		$CountMissingParamaters++;
		}
	}
	
#### Checking if all requested parameters were defined
	
	if ($CountMissingParamaters > 0) {
	print "\nThe following one-line-command parameters are needed, please provide them and try again:\n";
		foreach $parameter (sort keys %hashMissingParameters) {
		print "-$parameter\n";
		}
	die "\nThis program is exiting now by module LoadParameters::Parameters\n\n";
	}else{
	print "\nAll parameters needed were succesfully loaded\n";
	}

}

1;


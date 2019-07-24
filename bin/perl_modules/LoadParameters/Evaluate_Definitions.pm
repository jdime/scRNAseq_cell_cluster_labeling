
### Modules LoadParameters::Evaluate_Definitions and LoadParameters::Parameters
### work together to request and evaluate parameters for Perl programs

### They were inspired by libraries "argparse" from Python and "optparse" from R

### Module LoadParameters::Evaluate_Definitions evaluates that one-command-line parameters provided by the User
### meet the expected format (numeric, integer, strings, <comma> separated options, etc)

package LoadParameters::Evaluate_Definitions;
require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw(
Evaluate_parameter_values_series_of_digits
Evaluate_numb_cols
Evaluate_numb_rows
Evaluate_cutoffs_negative
Evaluate_cutoffs_positive
Evaluate_restrict_classes
Evaluate_value
Evaluate_file_exist
Evaluate_inflation
Evaluate_type_of_tips
Evaluate_horizontal_offset
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

######################################################
### If you ADD MORE subroutines to this module you'll need to pass input via LoadParameters::Parameters.pm
### Also make sure to add their main output (for external code) to @EXPORT
######################################################

######################################################
### This subroutine is for parameters whose format needs to be checked for specifit formats
######################################################

#### IMPORTANT: keep this row to START parsing of 'make_list_of_perl_dependencies_in_perl_programs.pl' script

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Evaluate_value {
my($parameterKey,$value) = @_;
	if ($parameterKey =~ /(^classes_type$)/) {
		unless ($value =~ /^(gmt|profile)$/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^cutoff_neg$)/) {
	### Negative values or NA or DEFAULT
		unless ($value =~ /(^-\d+$)|(^-\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^cutoff_ora$)/) {
	### Positive values or NA or DEFAULT
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^cutoff_pos$)/) {
	### Positive values or NA or DEFAULT
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^cutoff_print_p$)/) {
	### Positive values or NA or DEFAULT
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^cutoff_print_q$)/) {
	### Positive values or NA or DEFAULT
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^lwd$)/) {
	### Positive values or NA or DEFAULT
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)|(^NA$)|(^default$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^nperm$)/) {
	### Positive values. Don't allow 'NA'
		unless ($value =~ /(^\d+$)|(^\d+\.\d+$)/) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^p_correction_test$)/) {
		unless ($value =~ /(^BY$)|(^BH$)|(^bonferroni$)/) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^print_plot_ticks_labels_legend$)/) {
		unless ($value =~ /(^all$)|(^ticks$)|(^legend$)|(^label_tick$)|(^none$)|(^NA$)/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^use_graphics_device$)/) {
		unless ($value =~ /^(pdf|png|NA)$/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^use_ora$)/) {
		unless ($value =~ /(^p$)|(^pc$)/) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^use_universe$)/) {
		unless ($value =~ /(^u$)|(^i$)|(^n1$)|(^n2$)/) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}elsif ($parameterKey =~ /(^use_values_or_ranks$)/) {
		unless ($value =~ /^(values|ranks)$/i) {
		die "ERROR!!! Reading parameter '$parameterKey' value '$value' don't match expected format\n";
		}

	}else{
	die "ERROR!!! Reading parameter '$parameterKey' couldn't find evaluation parameters in LoadParameters::Evaluate_Definitions\n";
	}
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Evaluate_directory_exist {
my($parameterKey,$directory) = @_;
	if ($directory =~ /^none$/) {
	}elsif (-d $directory) {
	}elsif ($parameterKey =~ /(^path_outfiles)/) { ### these are directories for outputs, thus if don't exist will be created
	`mkdir -p $directory`;
	print "WARNING!!! created directory\n'$directory'\n";
	}elsif ($parameterKey =~ /(^hyperG_directory$)/) { ### this is an intermediate directory, thus just look that there is a defined string
		if ($directory =~ /\S/) {
		}else{
		die "ERROR!!! directory '$directory' was recognized as a string\n";
		}
	}else{
	die "ERROR!!! directory '$directory' was not found\n";
	}
}


##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Evaluate_file_exist {
my($parameterKey,$file) = @_;
	if ($file =~ /(^none$)|(^na$)/i) {
	}elsif (-f $file) {
	}else{
	die "ERROR!!! Reading parameter '$parameterKey' couldn't find file '$file'\n";
	}
}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Evaluate_were_defined {
my($parameterKey,$value) = @_;

	unless ($value =~ /\S/) {
	die "ERROR!!! Reading parameter '$parameterKey' string '$value' is not valid\n";
	}
}

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#### IMPORTANT: keep this row to END parsing of 'make_list_of_perl_dependencies_in_perl_programs.pl' script

1;


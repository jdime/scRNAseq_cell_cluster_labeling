#########################
### This module determines the HOME path and USER_NAME
###
### If running in a Linux the HOME will be '/home', in a Mac it will be '/Users'
###
### The resulting path to Users HOME is exported by $Users_home
### The resulting User name is exported by $DefaultUserName
#########################

package PathsDefinition::PathsToInputs;
require Exporter;
require AutoLoader;
@ISA = qw( Exporter AutoLoader );
@EXPORT = qw( $Users_home $DefaultUserName );

$User_home_path = `echo ~`;
chomp $User_home_path;
if ($User_home_path =~ /^(\/)(\S+)(\/)(\S+)/) {
$Users_home = $2;
$DefaultUserName = $4;
}

1;

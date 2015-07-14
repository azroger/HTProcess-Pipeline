#!/usr/bin/perl -w

use File::Slurp;
use warnings;

# written by Roger Barthelson
# This is used with HTProcess applications to embed png images directly 
# in a final html file, so that it can be opened directly in a browser
# from a data store file in the DE
# Dennis Roberts helped with some ideas for this

# embed png files in fastqc output, Roger Barthelson

$webfile = $ARGV[0];
$out_file = $ARGV[1];
if ( !defined $webfile || !defined $out_file  ) {
    die "Usage:  $0  webfile outfile \n";
}
open( WFIL, "$webfile" ) or die "Cannot open webfile\n";
open( OFIL, ">>$out_file" ) or die "Cannot open outfile\n";
while ( $lin = <WFIL> ) {
    chomp($lin);
    $lin =~ tr/\t/ /;
    @col = split( / /, $lin );
    my $col_count = scalar(@col);
    $max = $col_count - 1 ;
    for ( $i = 0 ; $i < $col_count ; $i++ ) {
    if ( $col[$i] =~ m/tick.png/ ) {
		my $translate_command = "openssl base64 -in Icons/tick.png -out tick.b64";
		system("$translate_command");
		$base64image = read_file('./tick.b64');
		print OFIL 'src="data:image/png;base64,';
		print OFIL "$base64image";
		print OFIL '"';
} elsif ( $col[$i] =~ m/warning.png/ ) {
		my $translate_command = "openssl base64 -in Icons/warning.png -out warning.b64";
		system("$translate_command");
		$base64image = read_file('./warning.b64');
		print OFIL 'src="data:image/png;base64,';
		print OFIL "$base64image";
		print OFIL '"';
} elsif ( $col[$i] =~ m/error.png/ ) {
	my $translate_command = "openssl base64 -in Icons/error.png -out error.b64";
	system("$translate_command");
	$base64image = read_file('./error.b64');
	print OFIL 'src="data:image/png;base64,';
	print OFIL "$base64image";
	print OFIL '"';
} elsif ( $col[$i] =~ m/fastqc_icon.png/ ) {
		my $translate_command = "openssl base64 -in Icons/fastqc_icon.png -out fastqc_icon.b64";
		system("$translate_command");
		$base64image = read_file('./fastqc_icon.b64');
		print OFIL 'src="data:image/png;base64,';
		print OFIL "$base64image";
		print OFIL '"';
} elsif ( $col[$i] =~ m/Images/ ) {
#$col[$i] =~ tr/src=//;
	@text = split (/=/, $col[$i] );
		my $translate_command = "openssl base64 -in $text[1] -out imagefile.b64";
		system("$translate_command");
		$base64image = read_file('./imagefile.b64');
		print OFIL 'src="data:image/png;base64,';
		print OFIL "$base64image";
		print OFIL '"';
} elsif ($i == $max) {
print OFIL "$col[$i]\n";
} else {
print OFIL "$col[$i] ";
}
}
}


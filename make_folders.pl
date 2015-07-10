#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my (@query_file1, @query_file2, @query_fileS, $null);
my $query_file1;
my $query_file2;
my $query_fileS;

GetOptions( "query_file1=s" => \@query_file1,
	    "query_file2=s" => \@query_file2,
	    "query_fileS=s" => \@query_fileS,
	    "output=s" => \$null );

#######################################
my $one = 'paired_reads1';
my $two = 'paired_reads2';
my $singles = 'single_reads';

if ( scalar(@query_file1) > 0 ) {
system( "mkdir $one" );
system( "mkdir $two" );
	for $query_file1 (@query_file1) {
	print "the query file is @query_file1 \n";    
 	system( "mv $query_file1 ./paired_reads1/" );   
    }

	for $query_file2 (@query_file2) {
	print "the query file is @query_file2 \n";
 	system( "mv $query_file2 ./paired_reads2/" );    
    }
}
if ( scalar(@query_fileS) > 0 ) {
system( "mkdir $singles" );
	for $query_fileS (@query_fileS) {
	print "the query file is @query_fileS \n";
 	system( "mv $query_fileS ./single_reads/" ); 
    }
}


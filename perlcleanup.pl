#!/usr/bin/perl

#use strict;
use warnings;

# used with htprocess_cuffdiff for doing some of the CummeRbund analsis
# written originally by Sheldon McKay (I think)
# requires also an R script, basic_plots.R

open( AFIL,  ">> CleanupOut" )  or die "Cannot open output file\n";

system "cuffdiff_get_goodstuff_test.pl 0.05";

open( IFIL,  "config.txt" );

	while ($lin = <IFIL>) {
	chomp($lin);

	@col = split (/\t/, $lin);
	push @labels,  "$col[0]" ;
#nrs++;
	}
	
#}
# Now plot a few graphs
my @pairs;
my %pair;
for my $i (@labels) {
    for my $j (@labels) {
        my $pair = join '.', sort ($i,$j);
        next if $pair{$pair}++ || $i eq $j;
        push @pairs, [$i,$j];
    }
}
my $r = <<END;
library(cummeRbund)
cuff <- readCufflinks()
png('../graphs/density_plot.png')
csDensity(genes(cuff))
dev.off()
END
;

for (@pairs) {
    my ($i,$j) = @$_;
    $r .= <<END;
png("../graphs/$i\_$j\_scatter_plot.png")
csScatter(genes(cuff),"$i","$j",smooth=T)
dev.off()
png("../graphs/$i\_$j\_volcano_plot.png")
csVolcano(genes(cuff),"$i","$j");
dev.off()
END
;

}

open RS, ">CUFFDIFFOUTPUT/basic_plots.R";
print RS $r;
close RS;
#close NAMES;
close AFIL;

system "mkdir graphs";
chdir "CUFFDIFFOUTPUT";
system "perl -i -pe 's/(\\S)\\s+?\\#/$1\\-/' *.*";
system "R --vanilla < basic_plots.R";

system "rm -f ../munged.gtf" if -e "../munged.gtf";


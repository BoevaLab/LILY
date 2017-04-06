#!/usr/bin/perl -w

#creates a .gff out of .bed caclulated with HMCan
#Copywrite Valentina Boeva, 2016
#This code cannot be made public


use strict;

my $SitesFilename ="";	
my $verbose = 0;
my $header = 0;
my $minscore = 0;
my $gffFile = 0;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f filename 	file with sites in BED format
    
    -----------------------------
    optional parameters:
    -v				for verbose
    -minScore			minimal score
    -head			if there is a header

    -o		filename	output file
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

my $ResFilename = "";

my $regFilename = "";

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }
    elsif ( $this_arg eq '-f') {$SitesFilename = shift @ARGV;}

    elsif ( $this_arg eq '-v') {$verbose = 1;}  
    elsif ( $this_arg eq '-head') {$header = 1;}

    elsif ( $this_arg eq '-o') {$gffFile = shift @ARGV;}
    elsif ( $this_arg eq '-minScore') {$minscore = shift @ARGV;}   
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ($gffFile eq "") {
    $gffFile = $SitesFilename.".gff";
}


open (FILE, "<$SitesFilename") or die "Cannot open file $SitesFilename to read!!!!: $!";

my %hash;
(<FILE>) if ($header);

my $count = 0;
my $countAccepted = 0;

while (<FILE>) {
	chomp;
	$count=$count+1;
	if (/(chr\S+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)/){
		if ($5>=$minscore) {
			$countAccepted=  $countAccepted+1;
			$hash{$1}->{$2}=$_;
		}
	}		
}
print "Total regions : $count\nAccepted:  $countAccepted\n" if $verbose;



close FILE;

open (OUT, ">$gffFile") or die "Cannot open file $gffFile to write!!!!: $!";

for my $chr (sort bychrom keys %hash) {
	for my $start (sort {$a <=> $b} keys %{$hash{$chr}}) {
		$_=$hash{$chr}->{$start};
		if (/(chr\S+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)/){
			print OUT join ("\t",$1,$4,"enhancer",$2,$3,$5,"+",".",$4), "\n";
		}
	} 
}


close OUT;
sub bychrom {
	my $c = $a;
	my $d = $b;
	$c=~s/chr//;
	$d=~s/chr//;
	$c=~s/X/23/;
	$d=~s/X/23/;
	$c=~s/Y/24/;
	$d=~s/Y/24/; 
	return $c <=> $d;
	
}

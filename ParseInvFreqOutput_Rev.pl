#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Statistics::Basic qw(:all);
use Number::Format;
my $en = new Number::Format(-decimal_point   => '.');

my %options=();
getopts("i:o:", \%options);
open(IN, "<", $options{i}) || die "no name file";
open(OUT, ">", $options{o}) || die "no input file";
my $gen=0;
my %hash;

while(<IN>)
{
	my @Read=split(/[\s\b\t\n]+/, $_);
	if ($_ =~ m/^Genera/)
	{
		for my $key (sort { $a <=> $b }keys %hash)
			{
			print OUT $gen, "\t",$hash{$key};
			}
		undef %hash;
		$gen=$Read[1];
	}
	else
	{
		$hash{$.}=$_;
	}
}	

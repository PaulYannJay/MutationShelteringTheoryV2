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

my %hash;

while(<IN>)
{
	my @Read=split(/[\s\b\t\n]+/, $_);
	if ($_ =~ m/^Genera/)
	{
		print OUT $Read[1], "\t";
		for my $key (sort { $a <=> $b }keys %hash)
			{
			my $MeanRec=mean(@{$hash{$key}});
			my $MeanFormat= $en->format_number($MeanRec);
			print OUT $MeanFormat, "\t";
			}
		print OUT "\n";
		undef %hash;
	}
	else
	{
	for my $i (1 .. $#Read)
		{	
		if (exists $hash{$i})
			{
			push @{$hash{$i}}, $Read[$i];
			}
		else
			{
			$hash{$i}=[ $Read[$i] ];
			}
		}
	}
}	

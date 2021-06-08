#!/usr/bin/perl

# Author: Beifang Niu

# This script gets pair-end reads recruitment information from
# FR-HIT results.

use strict;
use warnings;
use Carp;

my $end0      = shift;
my $end1      = shift;
my $endout0   = shift;
my $endout1   = shift;
my $insert    = shift;
my $deivation = shift;

unless ($end0) {
    print STDERR "\nUsage: frhit2pairend.pl <end0 fasta file> <end1 fasta file> <end0 FR-HIT result> <end1 FR-HIT result> [insert size:default 200] [deviation length: default 100]\n\n";
    print STDERR "Results file: FR-HIT.pair.out (pair-end recruitment information); FR-HIT.single.out (single-end recruitment information); \n\n";
    exit 0;
}

my (@nameend0, @nameend1, @outend0, @outend1, @conts);

# insert size
$insert = 200 if (($insert> 1000) || ($insert<0) || ($insert eq ""));
$deivation = 100 unless ($deivation);

my $maxinsert = $insert + $deivation;
my $mininsert = $insert - $deivation;

# expect valid insert size
$mininsert = 0 if ($mininsert < 0);

# load one end
open(IN,"$end0") || confess "could not open '$end0': $!";
my @entireFile = <IN>;
close(IN);

@nameend0 = map{ /^>(\S+)/ } @entireFile;

# load another end
open(IN,"$end1") || confess "could not open '$end1': $!";
@entireFile = <IN>;
close(IN);
@nameend1 = map{ /^>(\S+)/ } @entireFile;

# load output 1
open(IN,"$endout0") || confess "could not open '$endout0': $!";
@entireFile = <IN>;
close(IN);
@outend0 = map{ /(.*?)\n/ } @entireFile;

# load output 2
open(IN,"$endout1") || confess "could not open '$endout1': $!";
@entireFile = <IN>;
close(IN);
@outend1 = map{ /(.*?)\n/ } @entireFile;
my $pairs = scalar(@nameend0);

# output file
my $pairout="FR-HIT.pair.out";
my $singleout="FR-HIT.single.out";

open(PAIR, "> $pairout");
open(SINGLE, "> $singleout");

foreach my $i (0..$pairs)
{
    my $catpair="";
    my $catsingle="";
    my $na0=$nameend0[$i];
    my $na1=$nameend1[$i];

    while (1) 
    {
        if (scalar(@outend0)>0) 
        {
            my $ll=$outend0[0];
            my @t=split(/\t/,$ll);
            if ($na0 eq $t[0]) 
            {
                push(@conts,[$t[0],$t[8],$t[9],$ll,$t[6]]);
                shift(@outend0);
            }
            else { last;}
        }
        else { last; }
    }
    while (1) 
    {
        if (scalar(@outend1)>0) 
        {
            my $ll=$outend1[0];
            my @t=split(/\t/,$ll);
            if ($na1 eq $t[0]) 
            {
                push(@conts,[$t[0],$t[8],$t[9],$ll,$t[6]]);
                shift(@outend1);
            }
            else { last; }
        }
        else { last; }
    }
    @conts = sort{$a->[2] <=> $b->[2]} @conts;
    @conts = sort{$a->[1] cmp $b->[1]} @conts;

    while (1) 
    {
        if (scalar(@conts)>1) 
        {
            my $a=shift(@conts);
            my $b=$conts[0];

            if ($a->[0] eq $b->[0])
            {
                #print $a to single end result
                my $catsingle.=$a->[3]."\n";
            }
            elsif ($a->[4] eq $b->[4])
            {
                #print $a to sinlge end result
                my $catsingle.=$a->[3]."\n";
            }
            elsif ($a->[1] ne $b->[1])
            {
                #print $a to sinlge end result
                my $catsingle.=$a->[3]."\n";
            }
            else
            {
                my $span=abs($b->[2]-$a->[2]);
                if (($mininsert<=$span)&&($span<=$maxinsert)) 
                {
                    #print pair infor
                    my $catpair.=$a->[3]."\n";
                    $catpair.=$b->[3]."\n";
                    shift(@conts);
                }
                else
                {
                    #print $a to sinlge end result
                    my $catsingle.=$a->[3]."\n";
                }
            }
        }
        else
        {
            if (scalar(@conts)==1)
            {
                #output single end results
                my $a=shift(@conts);
                my $catsingle.=$a->[3]."\n";
            }

            last;
        }
    }

    @conts=();
    print PAIR $catpair;
    print SINGLE $catsingle;
}

close(PAIR);
close(SINGLE);


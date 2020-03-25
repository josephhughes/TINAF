#!/usr/bin perl

##################################################################### 
# A perlscript written by Joseph Hughes, University of Glasgow
#####################################################################
# Use this script to determine the depth of each amplicon from a BAM file 
# perl batch7c/Script/AmpliconDepthSystem.pl -bam batch7c/CVR100_bwa.primerclipped.bam -bedpe batch7c/primerV2.bedpe -out test

use strict;
use Getopt::Long; 
use warnings;
use POSIX;
#use Bio::SeqIO;
#use DateTime qw();
#my $date =DateTime->now->strftime('%Y%m%d');

my ($help,$bam,$out,$bedpe); 
&GetOptions(
      'bam:s'  => \$bam, # the bam file to check coverage for 
      'out:s' => \$out, # the output text tab delimited for each primer pair with coordinate start end and depth 
      'bedpe:s'  => \$bedpe, # the coordinates of each amplicon in bedpe format
      'help:s' => \$help,
           );
if ($help||!$bam||!$out||!$bedpe){
  print "Usage: perl batch7c/Script/AmpliconDepthSystem.pl -bam batch7c/CVR100_bwa.primerclipped.bam -bedpe batch7c/primerV2.bedpe -out test\n";
  print " -bam <bam> - the bam file to check\n";
  print " -out <txt> - output in text-tab delimited format\n";
  print " -bedpe <txt> - the coordinates of primers in bedpe format\n";
  print " -help        - Get this help\n";
  exit();
}

my %coord;
open(BED,"<$bedpe")||die "Can't open $bedpe\n";
while(<BED>){
  chomp($_);
  my @values=split(/\t/,$_);
  if ($values[0] eq $values[3]){
    $coord{$values[0]}{$values[6]}{"start"}=$values[1];
    $coord{$values[0]}{$values[6]}{"end"}=$values[5];
  }else{
    print "Different targets for $values[6]: $values[0] and $values[3]\n";
  }
}

open(OUT,">$out")||die "Can't open $out\n";
foreach my $target (keys %coord){
  foreach my $primer (sort {$coord{$target}{$a}{"start"} <=> $coord{$target}{$b}{"start"}} keys %{$coord{$target}}){
    my $cmd="samtools depth $bam -r $target:".$coord{$target}{$primer}{"start"}."-".$coord{$target}{$primer}{"end"};
    #print "$cmd\n";
    #my $depth=`$cmd`;
    open (my $pipe_fh, "$cmd |");
    my $sites=0;
    my $cumdepth=0;
    while(<$pipe_fh>){
      chomp($_);
      my @depth=split(/\t/,$_);
      $sites++;
      $cumdepth=$cumdepth+$depth[2];
    }
    print OUT "$primer\t".$coord{$target}{$primer}{"start"}."\t".$coord{$target}{$primer}{"end"}."\t".floor($cumdepth/$sites)."\n";
  }
}
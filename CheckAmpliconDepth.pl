#!/usr/bin perl

##################################################################### 
# A perlscript written by Joseph Hughes, University of Glasgow
#####################################################################
# Use this script to QC the depth of the amplicons

use strict;
use Getopt::Long; 
use warnings;
use POSIX;
#use Bio::SeqIO;
#use DateTime qw();
#my $date =DateTime->now->strftime('%Y%m%d');

my ($help,$ampdepth,$out,$depth); 
&GetOptions(
      'table:s'  => \$ampdepth, # the amplicon depth table for the run
      'out:s' => \$out, # the output text tab delimited for each primer pair with coordinate start end and depth 
      'd:i' => \$depth,
      'help:s' => \$help,
           );
if ($help||!$ampdepth||!$out){
  print "Usage: perl CheckAmpliconDepth.pl -table AmpliconDepth.txt -out test\n";
  print " -table <txt> - path to the amplicon depth table [compiled outputs from Sreenu's script]\n";
  print " -out <txt> - output in text-tab delimited format\n";
  print " -d <i> - minimum depth for amplicon\n";
  print " -help        - Get this help\n";
  exit();
}

open(TABLE,"<$ampdepth")||die "Can't open $ampdepth\n";
my $header=<TABLE>;
chomp($header);
my @colnames=split(/\t/,$header);
my %depth;
my %Nrange;
my %coord;
my %samples;
while(<TABLE>){
  chomp($_);
  my @values=split(/\t/,$_);
  my $neg_max_depth=0;
  for (my $i=3; $i<scalar(@values); $i++){
    if ($colnames[$i]=~/Neg|NC/ && $values[$i]>$neg_max_depth){
      $neg_max_depth=$values[$i];
    }elsif ($colnames[$i]!~/Neg|NC|Undetermined/){
      $depth{$values[0]}{$colnames[$i]}=$values[$i];
      $samples{$colnames[$i]}++;
    }
  }
  $depth{$values[0]}{"maxneg"}=$neg_max_depth;
  #print "Max negative depth for amplicon $values[0] is $neg_max_depth\n";
  $coord{$values[0]}{"label"}=$values[0];
  $coord{$values[0]}{"start"}=$values[1];
  $coord{$values[0]}{"stop"}=$values[2];
}

my %contamin;
my %failed;
for my $amplicon (keys %depth){
  foreach my $sample (keys %samples){
      if ($depth{$amplicon}{$sample}<(1.5*$depth{$amplicon}{"maxneg"}) && $depth{$amplicon}{$sample}>$depth){
        $contamin{$sample}{$amplicon}++;
#        print "$colnames[$i] $amplicon\n";
      }elsif ($depth{$amplicon}{$sample}<$depth){
        $failed{$sample}{$amplicon}++;
      }
   }
}
open(OUT,">$out")||die "Can't open $out\n";
foreach my $sample (keys %samples){
  print OUT "$sample\t";
  my @listfailed;
  my @listcontamin;
  my @Nregion;
  foreach my $failedamp (keys %{$failed{$sample}}){
    push(@listfailed,"Amp".$failedamp);
    #push(@Nregion,$Nrange{$failedamp});
  }
  foreach my $contaminamp (keys %{$contamin{$sample}}){
    push(@listcontamin,"Amp".$contaminamp);
    push(@Nregion,$coord{$contaminamp}{"start"}."-".$coord{$contaminamp}{"stop"});
  }
  
  
  
  print OUT join(",",@listfailed),"\t",join(",",@listcontamin),"\t",join(",",@Nregion),"\n";
}
#!/usr/bin perl

##################################################################### 
# A perlscript written by Joseph Hughes, University of Glasgow
#####################################################################
# perl CheckPrimerLocation.pl -fasta Wuhan_ref.fa.txt -primer nCoV-2019.csv > Primer_locations.txt

use strict;
use Getopt::Long; 
use warnings;
use Bio::SeqIO;
#use String::Approx 'amatch';

my ($help,$fasta,$csv,$tsv,$bedpe); 
&GetOptions(
	    'fasta:s'  => \$fasta, # the reference sequence
	    'csvprimer:s' => \$csv, # csv with the primer sequence in the 3rd column
	    'tsvprimer:s' => \$tsv, # tsv with the primer sequence in the 3rd column
	    'bedpe:s'  => \$bedpe, # bed file output
        "help:s" => \$help,
           );
           
if ($help||!$fasta||!$bedpe){
  print "Usage: perl CheckPrimerLocation.pl -fasta MN908947.3_ref.fat -csvprimer nCoV-2019.csv -bedpe primer.bedpe\n";
  print " -fasta <fasta> - the reference fasta sequence\n";
  print " -csvprimer <txt> OR -tsvprimer - list of primer sequences name in first column and sequence in third column\n";
  print " -bedpe <txt> - output coordinates of primers in bedpe format\n";
  print " -help        - Get this help\n";
  exit();
}



my %primers;
if ($csv){
  open(PRIMER,"<$csv") || die "Can't open $csv\n";

  while(<PRIMER>){
    chomp($_);
    my @cells=split(/,/,$_);
    $primers{$cells[0]}=$cells[2];
    my $revcomp = uc(reverse $cells[2]);
    $revcomp =~ tr/ATGC/TACG/;
    #print "$cells[2]\t$revcomp\n";
    $primers{$cells[0]."_revcomp"}=$revcomp;
  }
}

if ($tsv){
  open(PRIMER,"<$tsv") || die "Can't open $tsv\n";
  while(<PRIMER>){
    chomp($_);
    my @cells=split(/\t/,$_);
    $primers{$cells[0]}=$cells[2];
    my $revcomp = uc(reverse $cells[2]);
    $revcomp =~ tr/ATGC/TACG/;
    #print "$cells[2]\t$revcomp\n";
    $primers{$cells[0]."_revcomp"}=$revcomp;
  }
}

my $in  = Bio::SeqIO->new(-file => "$fasta" ,
                         -format => 'fasta');


my %primerpairs;
my $cnt=0;
while ( my $seq = $in->next_seq() ) {
  my $seq_str = uc($seq->seq());
  my $id = $seq->id();
  #print "$seq_str\n";
  for my $primer_name (keys %primers){
    if ($seq_str =~m/(.+)$primers{$primer_name}(.+)/){
      my $before=$1;
      my $after=$2;
      #print "$id has $primer_name $primers{$primer_name}\n";
      #print length($before),"\n";
      #print length($after),"\n";
      (my $number=$primer_name)=~s/nCoV-2019_(\d+)_.+/$1/;
      print "$number\t$primer_name\t$primers{$primer_name}\t",length($before)+1,"\t",length($before)+length($primers{$primer_name}),"\n";
      (my $pair=$primer_name)=~s/(.+)_(RIGHT|LEFT).*/$1/;
      (my $type=$primer_name)=~s/.+_(RIGHT|LEFT)(.*)/$1$2/;
      #print "$pair $type\n";
      if ($type=~/revcomp/){
        $primerpairs{$pair}{"reverse"}{"start"}=length($before)+1;
        $primerpairs{$pair}{"reverse"}{"end"}=length($before)+length($primers{$primer_name}); 
        $primerpairs{$pair}{"reverse"}{"id"}=$id;   
      }else{
        $primerpairs{$pair}{"forward"}{"start"}=length($before)+1;
        $primerpairs{$pair}{"forward"}{"end"}=length($before)+length($primers{$primer_name});  
        $primerpairs{$pair}{"forward"}{"id"}=$id;                
      }
      
    }
   # if (amatch($primers{$primer_name}, [ "2" ], $seq_str)){
   #   print "With up to 2 mismatches: $id has $primer_name $primers{$primer_name}\n";
   # }
  }
}
open(BED,">$bedpe")||die "Can't open $bedpe\n";
foreach my $pairs (keys %primerpairs){
  print BED "$primerpairs{$pairs}{'forward'}{'id'}\t$primerpairs{$pairs}{'forward'}{'start'}\t$primerpairs{$pairs}{'forward'}{'end'}\t";
  print BED "$primerpairs{$pairs}{'reverse'}{'id'}\t$primerpairs{$pairs}{'reverse'}{'start'}\t$primerpairs{$pairs}{'reverse'}{'end'}\t$pairs\n";
}

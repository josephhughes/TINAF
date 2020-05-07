#!/usr/bin/perl

# count the number of Ns in each sequence

# load module
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($help,$fasta);

GetOptions( "fasta=s" => \$fasta, # fasta sequence
		    "help" => \$help, );   # to get help

if(!($help)&&!($fasta))
 {
 print "Usage : CountNs.pl <list of arguments>, all arguments are necessary\n";
 print " -fasta <txt>     - a fasta file with 1 or more sequences\n";
 print " -help          - Get more detailed help\n";
 exit();
 }
 
if($help)
 {
 print "Cout the Ns in a sequence\n\n";

 print "---------------------------------------------------------------------\n";
 print "Usage : CountNs.pl <list of arguments>, all arguments are necessary\n";
 print " -fasta <txt>     - a fasta file with 1 or more sequences\n";
 print " -help          - Get more detailed help\n";

exit();
 }

my %seqs;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $fasta);
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  $seqs{$id}=$seq->seq();
}

foreach my $id (keys %seqs){
  my $count = () = $seqs{$id} =~ /\QN/gi;
  print "$id\t$count\n";
}


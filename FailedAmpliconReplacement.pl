#!/usr/bin/perl

# input is the alignment and the coordinates to replace according to a reference
# the reference name needs to be given
# either provide a coordinate on the command line which will replace it for all sequences (except ref)
# or tab-delimited file with the name of each sequence and the coordinate
# multiple coordinates acan be specified by start1-end1,start2-end2
# the output is the alignment with the Ns in the correct positions for each of the sequences in the text-tab or all if specified on the command line

# load module
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($aln, $help, $alnout, $coord,$ref,$fasta);

GetOptions( "aln=s" => \$aln, # alignment file
	    	"ref=s" => \$ref, # the reference sequence name
	    	"coord=s" => \$coord, # single coordinate
	    	"alnout=s" => \$alnout, # alignment output 
	    	"fasta=s"  => \$fasta, # fasta sequence output without the gaps
		    "help" => \$help, );   # to get help

if(!($help)&&!($coord)||!($aln)||!($alnout)||!($ref))
 {
 print "Usage : FailedAmpliconReplacement.pl <list of arguments>, all arguments are necessary\n";
 print " -aln <txt>     - an alignment with IDs that match up with the metadata\n";
 print " -ref <txt>     - name of the reference sequence the coordinates are from\n";
 print " -coord <txt>   - coordinates in the format start1-end1,start2-end2, all sequences in the alignment will be replaced with these coordinates\n";
 print "                OR a text-tab file with the name of the sequence in the first column and the coordinates to replace in the second column\n";
 print " -alnout <txt>  - the output alignment with the Ns in failed regions\n";
 print " -fasta <txt>- the output in fasta unaligned\n";
 print " -help          - Get more detailed help\n";
 exit();
 }
 
if($help)
 {
 print "To replace defined regions on sequences with Ns for failed amplicons\n\n";

 print "---------------------------------------------------------------------\n";
 print "Usage : FailedAmpliconReplacement.pl <list of arguments>, all arguments are necessary\n";
 print " -aln <txt> - an alignment with IDs that match up with the metadata\n";
 print " -ref <txt>    - name of the reference sequence the coordinates are from\n";
 print " -coord <txt> - coordinates in the format start1-end1,start2-end2, all sequences in the alignment will be replaced with these coordinates\n";
 print " OR a text-tab file with the name of the sequence in the first column and the coordinates to replace in the second column\n";
 print " -alnout <txt> - the output alignment with the Ns in failed regions\n";
 print " -fasta <txt>- the output in fasta unaligned\n";
 print " -help          - Get more detailed help\n";

exit();
 }

my %seqs;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $aln);
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  $seqs{$id}=$seq->seq();
}

if (!$seqs{$ref}){
  print "$ref reference identifier is not in the alignment\n";
  exit();
}

if (-e $coord){
  my (%coord);
  # print "$coord is a file\n";
  open(COORD,"<$coord")||die "Can't open $coord\n";
  while(<COORD>){
    chomp($_);
    my @values=split(/\t/,$_);
    $coord{$values[0]}=$values[1];
    #print "$values[0]\n";
  }
  close(COORD);
  foreach my $id (keys %coord){
    my @regions=split(/,/, $coord{$id});
    #print "@regions\n";
    my @seq=split(//,$seqs{$id});
    my @refseq=split(//,$seqs{$ref});

    foreach my $region (@regions){
      my ($start,$end)=split(/-/,$region);
      my $refgap=0;
      for (my $i=0; $i<scalar(@seq); $i++){
        if ($refseq[$i] eq "-"){
          $refgap++;
        }
        my $position=$i+1-$refgap;
        if ($position<=$end && $position>=$start){
          if ($seq[$i] ne "-"){
            substr $seqs{$id}, $i, 1, "N";
          }
        }
      }
    }
    #print "$id\n$seqs{$id}\n";
  }
  
}else{#  in this case all the sequences get changed
  my @regions=split(/,/,$coord);
  foreach my $id (keys %seqs){
    my @seq=split(//,$seqs{$id});
    my @refseq=split(//,$seqs{$ref});
    foreach my $region (@regions){
      my ($start,$end)=split(/-/,$region);
      my $refgap=0;
      for (my $i=0; $i<scalar(@seq); $i++){
        if ($refseq[$i] eq "-"){
          $refgap++;
        }
        my $position=$i+1-$refgap;
        if ($position<=$end && $position>=$start){
          if ($seq[$i] ne "-"){
            substr $seqs{$id}, $i, 1, "N";
          }
        }
      }
    }
    #print "$id\n$seqs{$id}\n";
  }
}

foreach my $id (keys %seqs){
  my $count = () = $seqs{$id} =~ /\QN/gi;
  print "$id\t$count\n";
}

if ($alnout){
  open(OUTALN,">$alnout")||die "Can't open $alnout\n";
  foreach my $id (keys %seqs){
    print OUTALN ">$id\n$seqs{$id}\n";
  }
}
if ($fasta){ # strip all gaps
  open(FASTA,">$fasta")||die "Can't open $alnout\n";
  foreach my $id (keys %seqs){
    my $sequence=$seqs{$id};
    $sequence=~s/\-//g;
    print FASTA ">$id\n$sequence\n";
  }

}

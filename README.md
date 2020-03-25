# TINAF
Scripts for processing metadata from GISAID and reconstructing phylogeny

## Illumina read processing

### Processing of fastq files

Bash script for looping through a directory of fastq files and creating a bam file. Created to process amplicon sequencing.

```
loopilluminaamplicon.sh /path/to/fastq_files/ ref.fa primer.bedpe & 
```

Inputs provided in the following order:
1. The directory to the fastq files (the script is expecting files formatted as R1_001.fastq.gz and R2_001.fastq.gz)
2. The reference fasta file
3. A list of primers in BEDPE format used for trimming 

The script goes through the following steps:
1. Trimming of Illumina adapters and porr quality reads using TrimmomaticPE (unpaired reads are ignored for later steps)
2. Mapping to the references using bwa and producing a sorted bam file
3. Indexing the bam file
4. Trimming of primers using bamclipper.sh
5. Summary statistics using weeSAM
6. Removal of unmapped reads
7. Indexing the bam file
8. Consensus calling with a minimum coverage of 5

### Amplicon depth

Perl script for calculating the average depth for each of the amplicons.

```
perl AmpliconDepthSystem.pl -bam input.bam -bedpe primer.bedpe -out depth.txt
```

Input:
1. The bam file
2. the coordinates of the primer pairs in BEDPE format

Output:
A text-tab file with the name of the primer, the amplicon start and end and the average depth of that amplicon.

# This code can run FASTQC on multiple fastq files presented as a list
# Author Candida Vaz
#!/usr/bin/perl
use strict;
use warnings;

my$list=$ARGV[0];
open(LIN,"$list");
my@larr=<LIN>;
chomp(@larr);

foreach my$file(@larr)
 {
   print "running fastqc on $file\n";
   system("fastqc $file");   
 }

close (LIN);

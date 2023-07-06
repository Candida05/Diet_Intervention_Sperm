# This code runs Trimmomatic on fastq files
# Author Candida Vaz
#!/usr/bin/perl
use strict;
use warnings;

my$list=$ARGV[0];#list of the raw files to run Trimmomatic
open(IN,"$list");
my@arr=<IN>;
chomp(@arr);
my$file=();
my$name=();
my$ext=();

foreach $file(@arr)
  {
     (my$name,my$ext)=split(/\./,$file);
     print "Running Trimmomatic for $name.$ext\n";
     system ("java -jar /home/candidavaz/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar SE $file $name-trim-orig.fastq ILLUMINACLIP:/home/candidavaz/Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE-mod.fa:2:30:10 LEADING:3 TRAILING:3 CROP:46 SLIDINGWINDOW:4:15 MINLEN:10\n"); 
     print "done trimming $file\n"; 
     print "*******************\n";   
  }

close(IN);

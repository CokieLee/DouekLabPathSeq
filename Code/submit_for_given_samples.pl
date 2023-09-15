#!/usr/bin/perl
#
#
use strict;

if($#ARGV < 6){die" Usage: Code_Path InputDataPath SampleName Start_number End_number (DNA/RNA) OutputPath\n"}

my($codePath,$inputPath,$sampName,$start,$end,$type,$outputPath)=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6]);
print "running PathSeq for $end $start..$end samples\n";


open(LIST,"$inputPath/file_list")|| die "could not open list\n";

my @list = <LIST>;
 if($start>$#list){print "$start and $end are more than the number of samples\n"; exit;}
 if($end > $#list){print "running it till $end";$end=$#list}
 

print "file pairs $#list\n";

my @a=($start..$end);

for my $i (@a)
{
  chomp($list[$i]);
  print "$list[$i]\n";
  my$left = $list[$i];
  my$right= $left;
  $right=~ s/_R1_/_R2_/;
  print "$left,$right,$list[$i]\n";

# print "$i$left Read $right Read files are submitted\n";
  my $run_cmd = "qsub $codePath/PathSeqSubmitter.sh $codePath $sampName $inputPath/$left $inputPath/$right $type 300 $outputPath";
  print "$run_cmd\n";
 system("$run_cmd");
 }
close LIST;
print "$#list $list[$#list]\n";
#print @list;

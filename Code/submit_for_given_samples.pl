#!/usr/bin/perl
#
#
use strict;

if($#ARGV < 2){die" Usage: Start_number End_number (DNA/RNA)\n"}


my($start,$end,$type)=($ARGV[0],$ARGV[1],$ARGV[2]);
print "running PathSeq for $end $start..$end samples\n";


open(LIST,"../Input_Data/file_list")|| die "could not open list\n";

my @list = <LIST>;
 if($start>=$#list){print "$start and $end are more than the number of samples\n"; exit;}
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
  my $run_cmd = "qsub PathSeqSubmitter.sh 2023020_EVD68_NS Input_Data/$left Input_Data/$right $type 300";
  print "$run_cmd\n";
 system("$run_cmd");
 }
close LIST;
print "$#list $list[$#list]\n";
#print @list;

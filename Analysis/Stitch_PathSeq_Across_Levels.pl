#!/usr/bin/perl
#
# This script takes in the merged taxonomic and quantification information
# for each sample of a batch of data analyzed with CZID. It combines
# the data across samples by taxonomic level, and reouputs a file for each taxonomic level.
# Input: run script in a directory of all the tpm and pseudocounts files
# output for all samples in a pathseq analysis run.
#
# TODO: modify this script so it can be given an input and output directory

$pwd=`pwd`;
print "$pwd";

@types = qw(class family genus kingdom order phylum species);
# @types = qw(species);


@samples = `ls -1 $pwd`;
chomp(@samples);

foreach $type(@types)
{
  my (%tpms, %pseudocounts, @sample_list);

  foreach $sample(@samples)
  {
    ## Store pseudocounts and tpm data for each sample

    if($sample =~ /Sample_/)
    {
      push(@sample_list,$sample);
      $string.=",$sample";
      $sample_path = "$sample/RNA_merge_TaxAndQuant";
      my $baseName = $sample;
      $baseName =~ s/Sample_//;
      if(!-e "$sample_path"){print "$sample_path files don't exist\n"}
      if(-e "$sample_path")
      {
        $file1= `ls -1 $sample_path/RNA_$type\_$baseName\_tpm.csv`;
        $file2= `ls -1 $sample_path/RNA_$type\_$baseName\_pseudocounts.csv`;
        open(FILE1,$file1);
        while(<FILE1>)
        {
          chomp($_);
          ($descript,$number)=split(/\,/,$_);
          #print "$number\t$descript\n";
          $tpms{$descript}{$sample}=$number;
        }
        close FILE1;
        open(FILE2,$file2);
        while(<FILE2>)
        {
          chomp($_);
          ($descript,$number)=split(/\,/,$_);
          #print "$number\t$descript\n";
          $pseudocounts{$descript}{$sample}=$number;
        }
        close FILE2;
      }
    }
  }
    ## output stored pseudocounts and tpm data to files
    open(TYPE1,">RNA_$type\_tpm.csv");
    print TYPE1 "$type$string\n";
    foreach $key (keys %tpms)
    {
          print TYPE1 $key;
          foreach $sample(@sample_list)
            {
              $tpm=$tpms{$key}{$sample};
              if(!$tpms{$key}{$sample}){$tpm=0}
              print TYPE1 ",$tpm";
            }
          print TYPE1 "\n";
    }

    close TYPE1;
    open(TYPE2,">RNA_$type\_pseudocounts.csv");
    print TYPE2 "$type$string\n";
    foreach $key (keys %pseudocounts)
    {
          print TYPE2 $key;
          foreach $sample(@sample_list)
            {
              $psc=$pseudocounts{$key}{$sample};
              if(!$pseudocounts{$key}{$sample}){$psc=0}
              print TYPE2 ",$psc";
            }
          print TYPE2 "\n";
    }
    close TYPE2;

  undef %tpms;
  undef %pseudocounts;
  #   $files = system("ls -1 \"$path/$j/DNA_salmon_quant\"");
  #   print $files;
  $string = "";
}
print "@types\n";
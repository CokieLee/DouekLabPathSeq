#!/usr/bin/perl
#
#
$pwd=`pwd`;
print "$pwd";

@types = qw(class family genus kingdom order phylum species);
print "Stitching files across @types for each timepoint\n";

$path = $pwd;
$path =~ s/scripts//;


print "$path";
chomp($path);
open(LIST,"$path/list_of_conditions")||{die "could not open $path/list_of_conditions file\n"};

while(<LIST>)
{
   $subdir=$_;
   $pathname=$path."$subdir";
   print"$subdir";
   @list = `ls -1 $pathname`;
  $count=0;
  foreach $type(@types)
      {

   foreach $file(@list)
       {
          if($file =~ /Sample_/)
             {
               $count++;
               $index{$file}=$count;
                print "$count\t$file";
                 $quantpath= "$pathname/$file/RNA_salmon_quant";
                 if(!-e "$quantpath"){print "$quantpath does not exist"}
                 if(-e "$quantpath")
                   {
                       $file1= `ls -1 $quantpath/RNA_$type\_tpm.csv`;
                       $file2= `ls -1 $quantpath/RNA_$type\_pseudocounts.csv`;
                        
                       open(FILE1,$file1); 
                       while(<FILE1>)
                            {
                                chomp($_);
                                ($descript,$number)=split(/\,/,$_);
                                #print "$number\t$descript\n";
                                 $tpms{$descript}{$count}=$number;
                             
                             }
        
                          close FILE1;  
                       open(FILE2,$file2); 
                       while(<FILE2>)
                            {
                                chomp($_);
                                ($descript,$number)=split(/\,/,$_);
                                #print "$number\t$descript\n";
                                $pseudocounts{$descript}{$count}=$number;
                             
                            }
                    close FILE2;  
                    open(TYPE1,">$subdir_$type\_tpm.csv");
                    print TYPE1 "$type$string\n";
             foreach $key (keys %tpms)
              {
                  print TYPE1 $key;
                  foreach $k(77..150)
                   {
                       $tpm=$tpms{$key}{$k};
                       if(!$tpms{$key}{$k}){$tpm=0}
                       print TYPE1 ",$tpm";
                  }
                 print TYPE1 "\n";
               }
             close TYPE1;
             open(TYPE2,">DNA_$type\_pseudocounts.csv");
             print TYPE2 "$type$string\n";
             foreach $key (keys %pseudocounts)
              {
                  print TYPE2 $key;
                  foreach $l(77..150)
                   {
                       $psc=$pseudocounts{$key}{$l};
                       if(!$pseudocounts{$key}{$l}){$psc=0}
                       print TYPE2 ",$psc";
                  }
                 print TYPE2 "\n";
               }
             close TYPE2;
          } 

                 } 
             }
        }
}
close LIST;


foreach $type(@types)
         {
           foreach $j(77..150)
              {
                  $string.=",DNA_$j";
                  $sample_path = "$path/$j/DNA_salmon_quant";
                  if(!-e "$sample_path"){print "$j files don't exist\n"}
                  if(-e "$sample_path"){
                    $file1= `ls -1 $sample_path/DNA_$type\_tpm.csv`;
                    $file2= `ls -1 $sample_path/DNA_$type\_pseudocounts.csv`;
                    open(FILE1,$file1); 
                    while(<FILE1>)
                        {
                              chomp($_);
                             ($descript,$number)=split(/\,/,$_);
                             #print "$number\t$descript\n";
                             $tpms{$descript}{$j}=$number;
                             
                        }
        
                    close FILE1;  
                    open(FILE2,$file2); 
                    while(<FILE2>)
                        {
                              chomp($_);
                             ($descript,$number)=split(/\,/,$_);
                             #print "$number\t$descript\n";
                             $pseudocounts{$descript}{$j}=$number;
                             
                        }
        
                    close FILE2;  
             } 
       open(TYPE1,">DNA_$type\_tpm.csv");
       print TYPE1 "$type$string\n";
       foreach $key (keys %tpms)
        {
            print TYPE1 $key;
            foreach $k(77..150)
             {
                 $tpm=$tpms{$key}{$k};
                 if(!$tpms{$key}{$k}){$tpm=0}
                 print TYPE1 ",$tpm";
            }
           print TYPE1 "\n";
         }
       close TYPE1;
       open(TYPE2,">DNA_$type\_pseudocounts.csv");
       print TYPE2 "$type$string\n";
       foreach $key (keys %pseudocounts)
        {
            print TYPE2 $key;
            foreach $l(77..150)
             {
                 $psc=$pseudocounts{$key}{$l};
                 if(!$pseudocounts{$key}{$l}){$psc=0}
                 print TYPE2 ",$psc";
            }
           print TYPE2 "\n";
         }
       close TYPE2;
    } 

       undef %tpms;
       undef %pseudocounts;
 #   $files = system("ls -1 \"$path/$j/DNA_salmon_quant\"");
 #   print $files;
 $string = "";
}
print "@types\n";


#!/usr/bin/perl

sub trim($);

   $c_f_dir="./";
   $l_c_f_file=$ARGV[0];
   $l_cf_file=$c_f_dir.$l_c_f_file;

   @lig=();
   $l_counter=0;

   open (IN, $l_cf_file);
   while($line=<IN>)
   {
      chomp $line;

      $chain=substr($line,21,1);
      $chain=trim($chain);

      $residue_type=trim(substr($line,17,3));
      $residue_num=trim(substr($line,22,4));
      $contact_f=trim(substr($line,62,4));

      $key=$residue_num.".".$residue_type.".".$chain;

      $l_checker=0;
      for($l=0;$l<$l_counter;$l++)
      {
         if($key eq $lig[$l][0])
         {
            if($contact_f > $lig[$l][1])
            {
               $lig[$l][1]=$contact_f;
            }
            $l_checker=1;
            $lig[$l][2]=$lig[$l][2]+$contact_f;
         }
      }
      if($l_checker==0)
      {
         $lig[$l_counter][0]=$key;
         $lig[$l_counter][1]=$contact_f;
         $lig[$l_counter][2]=$contact_f;
         $l_counter++;
      }
   
   }
   close IN;

$lig_filename=$ARGV[0].".RCF";

open(OUT1, ">$lig_filename");

print (OUT1 "Residue_number\tResidue_type\tChain\tRCF\n");

for($g=0;$g<$l_counter;$g++)
{
   @sample=();
   @sample=split(/\./,$lig[$g][0]);
   print (OUT1 "$sample[0]\t$sample[1]\t$sample[2]\t$lig[$g][2]\n");
}
close OUT1;

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

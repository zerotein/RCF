#!usr/bin/perl

if(@ARGV!=2){
        print "\nsyntax:>perl RCF.pl zdock_out_file rotational_sampling_mode\n";
        print "\nFor rotational_sampling_mode type FG for 54K fine sampling mode\n";
        print "Or type CG for 3600 coarse sampling mode\n";
        print "\nExample,\n";
        print ">perl RCF.pl 1A2K.zd3.0.cg.out CG\n";
        exit;
        }

$angle="";

if($ARGV[1] eq "CG")
{
   $angle="euler.15";
}
elsif($ARGV[1] eq "FG")
{
   $angle="euler.6";
}
else
{
   print "\nPlease choose the either FG or CG for the rotational sampling mode\n";
   print "\nExample, RCF.pl 1A2K.zd3.0.cg.out CG\n";
   exit;
}

$r="";
$l="";

open(IN, $ARGV[0]);
$counter=0;
while($line=<IN>)
{
   chomp $line;
   @sample=();
   if($counter==2)
   {
      @sample=split(/\t/,$line);
      $r=$sample[0].".$ARGV[1]";
   }
   elsif($counter==3)
   {      
      @sample=split(/\t/,$line);
      $l=$sample[0].".$ARGV[1]";
   }
   $counter++;

   if($counter > 5)
   {
      break;
   }
}
close IN;

$command="source/reformat $ARGV[0] source/$angle 3.0";
system($command);

$hout=$ARGV[0].".hout";

$command="source/contact_frequency -i $hout -o $ARGV[1]";
system($command);

$command="perl source/Cal_RCF.pl $r; perl source/Cal_RCF.pl $l";
system($command);

$command="rm $hout $r $l";
system($command);

print "$ARGV[0] is done!!\n";

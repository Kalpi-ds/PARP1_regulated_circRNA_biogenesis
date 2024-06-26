#!/usr/bin/perl -w
use strict;

# SAMPLE1: PJ
# SAMPLE2: WT

my $inFN = "circRNAs.pVal.FDR.txt";
open(INFILE, $inFN) || die;
open(BOTH, ">circRNAs.pVal.FDR.txt.BOTH") || die;
open(WT, ">circRNAs.pVal.FDR.txt.SAMPLE2") || die;
open(KD, ">circRNAs.pVal.FDR.txt.SAMPLE1") || die;
open(WTONLY, ">circRNAs.pVal.FDR.txt.SAMPLE2ONLY") || die;
open(KDONLY, ">circRNAs.pVal.FDR.txt.SAMPLE1ONLY") || die;
my $hdr = <INFILE>;
print BOTH $hdr;
print WT $hdr;
print KD $hdr;
print WTONLY $hdr;
print KDONLY $hdr;
chomp($hdr);
 my @wds = split(/\t/, $hdr);
 my $numWds = @wds;
 for(my $i = 0; $i < $numWds; $i++) { 
    print "$i\t$wds[$i]\n";
}
print "KD\tWT\n";
while(defined(my $line = <INFILE>)) { 
   my @wds = split(/\t/, $line);
   print "$wds[12]\t$wds[14]\n";
   my @wds1 = split(/\,/, $wds[12]);
   my $numWds1 = @wds1;
   my $sumKD = 0;
   for(my $i = 0; $i < $numWds1; $i++) { 
      $sumKD += $wds1[$i];
   }

   @wds1 = split(/\,/, $wds[14]);
   $numWds1 = @wds1;
   my $sumWT = 0;
   for(my $i = 0; $i < $numWds1; $i++) { 
      $sumWT += $wds1[$i];
   }
   if($sumWT > 0) { 
      print WT "$line";
      if($sumKD == 0) { 
         print WTONLY "$line";
      }
   }
   if($sumKD > 0) { 
      print KD "$line";
      if($sumWT == 0) { 
         print KDONLY "$line";
      }
   }
   if(($sumWT > 0) && ($sumKD > 0)) {
      print BOTH "$line";
   }
}
close(INFILE);
close(BOTH);
close(WT);
close(KDONLY);
close(WTONLY);
close(KD);

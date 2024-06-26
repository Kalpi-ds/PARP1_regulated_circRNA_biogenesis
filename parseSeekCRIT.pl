#!/usr/bin/perl -w
use strict;

## Author: Eric Rouchka
## Date:   4-14-2022

my %genomeHASH = getGenomeHash("/bio/data/Genomes/Dm/BDGP6/BDGP6.fa");
my $seekCRITDIR = "/bio/home/kalina/KBRIN03XX-Yvonne-BECCA/SeekCRIT";
my @comparisonARR = ("DKDvsDWT", "DPJvsDWT", "OKDvsOWT_4");
my $numComparisons = @comparisonARR;
my %typeHASH;

for(my $i = 0; $i < $numComparisons; $i++) {
   my $currComp = $comparisonARR[$i];
   my $currCompFN = $seekCRITDIR . "/seekCRIT_out_" . $currComp . "/seekCRIT_output/circRNAs.pVal.FDR.txt";
   open(INFILE, $currCompFN) || die("Error opening $currCompFN for reading");
   my $intronFN = "$currComp.INTRONIC.fasta";
   my $exonFN   = "$currComp.EXONIC.fasta";
   my $allFN    = "$currComp.ALL.fasta";
   open(INTRONFN, ">$intronFN") || die;
   open(EXONFN, ">$exonFN") || die;
   open(ALLFN, ">$allFN") || die;

   my $hdr = <INFILE>;
   chomp($hdr);
   my @wds = split(/\t/, $hdr);
   my $numWds = @wds;
   my $exonOnly = 0;
   my $intronOnly = 0;
   my $intronExon = 0;

   for(my $j = 0; $j < $numWds; $j++) { 
#      print "$j\t$wds[$j]\n";
   }
   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      my @wds = split(/\t/, $line);
      my $chr = $wds[0];
      my $circStart = $wds[1];
      my $circEnd   = $wds[2];
      my $strand    = $wds[3];
      my $exonSizes = $wds[5];
      my $exonOffsets = $wds[6];
      my $circType = $wds[7];
      $typeHASH{$circType} = $circType;
      if($strand eq "-") { 
         $circStart += 1; ## ERROR IN OFFSET OF seekCRIT
      }
      my @intronBeg;
      my @intronEnd;
      my @exonBeg;
      my @exonEnd;
      my @sizeARR = split(/\,/, $exonSizes);
      my @offsARR = split(/\,/, $exonOffsets);
      my $numExons = @sizeARR;  ## Note there is an extra , at the end
      for(my $j = 0; $j < $numExons; $j++) { 
          my $currB = $offsARR[$j];
          my $currE = $currB + $sizeARR[$j] - 1;
          push(@exonBeg, $currB);
          push(@exonEnd, $currE);
      }
      my $sLen = $circEnd - $circStart + 1;
      if($exonBeg[0] > 0) { 
         push(@intronBeg, 0);
         push(@intronEnd, $exonBeg[0] - 1);
      }
      for(my $j = 1; $j < $numExons; $j++) { 
         push(@intronBeg, $exonEnd[$j - 1] + 1);
         push(@intronEnd, $exonBeg[$j] - 1);
      }
      my $lastIntronBeg = $exonEnd[$numExons - 1] + 1;
      my $lastIntronEnd = $sLen;
      if($lastIntronEnd > $lastIntronBeg) { 
         push(@intronBeg, $lastIntronBeg);
         push(@intronEnd, $lastIntronEnd);
      }
      my $intronSizes;
      my $intronOffsets;
      my $numIntrons = @intronBeg;
      for(my $j = 0; $j < $numIntrons; $j++) { 
          my $currSize = $intronEnd[$j] - $intronBeg[$j] + 1;
          $intronSizes .= "$currSize,";
          $intronOffsets .= "$intronBeg[$j],";
      }    

      my $chrSeq = $genomeHASH{$chr};
      my $circSeqALL = getSubSeq($chrSeq, $circStart, $circEnd);
      my $exonSeq = getExonSeq($circSeqALL, $exonSizes, $exonOffsets);
      my $intronSeq = "";
      if($numIntrons > 0) { 
         $intronSeq = getExonSeq($circSeqALL, $intronSizes, $intronOffsets);
      }

      if($strand eq "-") { 
         $circSeqALL = reverseComplement($circSeqALL); 
         $exonSeq = reverseComplement($exonSeq);
         $intronSeq = reverseComplement($intronSeq);
      }
      if($circType eq "ciRNA") { ## circular intronic
         print INTRONFN ">$line\n$circSeqALL\n";
         $intronOnly++;
      }
      else {
         if(length($exonSeq) > 0) { 
            print EXONFN ">$line\n$exonSeq\n";
            if(length($intronSeq) == 0) { 
               $exonOnly++;
            }
            else {
               $intronExon++;
            }
         }
         if(length($intronSeq) > 0) { 
            print INTRONFN ">$line\n$intronSeq\n";
         }
      }
      print ALLFN ">$line\n$circSeqALL\n";
     
#      print "$line\n";
#      print "ALLSEQ: $circSeqALL\nEXON: $exonSeq\nINTRON: $intronSeq\n";
   }
#0       chrom
#1       circRNA_start
#2       circRNA_end
#3       strand
#4       exonCount
#5       exonSizes
#6       exonOffsets
#7       circType
#8       geneName
#9       isoformName
#10      exonIndexOrIntronIndex
#11      FlankingIntrons
#12      CircularJunctionCount_Sample_1
#13      LinearJunctionCount_Sample_1
#14      CircularJunctionCount_Sample_2
#15      LinearJunctionCount_Sample_2
#16      PBI_Sample_1
#17      PBI_Sample_2
#18      deltaPBI(PBI_1-PBI_2)
#19      pValue
#20      FDR
   print "EXON: $exonOnly; INTRON: $intronOnly; BOTH: $intronExon\n";

   close(INFILE);
   close(ALLFN);
   close(EXONFN);
   close(INTRONFN);
}

#--------------------------------------------------------------------------------------
sub getExonSeq {
   my($s, $siz, $offs) = @_;

   my $exonSeq = "";   
   my @sizeARR = split(/\,/, $siz);
   my @offsARR = split(/\,/, $offs);
   my $numExons = @sizeARR;  ## Note there is an extra , at the end
   for(my $i = 0; $i < $numExons; $i++) { 
      my $currS = substr($s, $offsARR[$i], $sizeARR[$i]);
      $exonSeq .= $currS;
   }
   return($exonSeq);
} 
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
sub getSubSeq {
   my($s, $beg, $end, $strand) = @_;

   my $seqLen = $end - $beg + 1;
   my $subSeq = substr($s, $beg - 1, $seqLen);   ## subtract 1 for 1-based array
   return($subSeq);
}
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
sub reverseComplement {
   my($s) = @_;

   $s = uc($s);  ## Convert to uppercase
   $s =~ s/A/t/g;
   $s =~ s/C/g/g;
   $s =~ s/G/c/g;
   $s =~ s/T/a/g;
   $s = reverse($s);
   $s = uc($s);
   return($s);
} 
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
sub getGenomeHash {
   my($inFN) = @_;
   open(INFILE, $inFN) || die("Error opening $inFN for reading\n");

   my %h;
   my $seq = "";
   my $ID = "";

   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      if($line =~ /^\>/) {
         if($ID ne "") { 
            $h{$ID} = $seq;
         }
         $ID = $line;
         $ID =~ s/^\>//g;
         $seq = "";
      }
      else {
         $seq .= $line;
      }
   }
   close(INFILE);
   return(%h);
}
#--------------------------------------------------------------------------------------


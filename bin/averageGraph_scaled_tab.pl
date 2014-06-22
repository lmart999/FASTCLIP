#!/usr/bin/perl
#
###########################################################
#
# This script is used to generate data for average graph.
#
###########################################################



# step 1, readin coordinate info
# and determine which genes to include

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;

if (!@ARGV) {
  print "command line: perl averageGraph_scaled_tab.pl /seq/chromosome/mm9/mm9_refseq_genes_UTR_annotation.txt treat_tabfile control_tabfile output_filename\n";
  exit;
}

$annofile = $ARGV[0];
$treatfile = $ARGV[1];
$controlfile = $ARGV[2];
$filename = $ARGV[3];



$thrd_abundance = 0;
$thrd_top = 100000;
$noZero = 0;
$cut = 0;


$utr5_scale = 200;
$cds_scale = 200;
$utr3_scale = 200;
$total_scale = $utr5_scale + $cds_scale + $utr3_scale;

my %utr5_start;
my %utr5_end;
my %utr3_start;
my %utr3_end;
my %cds_start;
my %cds_end;

open (in, "<$annofile");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  if ($data[1] eq "5UTR") {
    $utr5_start{$data[0]} = $data[2];
    $utr5_end{$data[0]} = $data[3];
  }
  if ($data[1] eq "3UTR") {
    $utr3_start{$data[0]} = $data[2];
    $utr3_end{$data[0]} = $data[3];
  }
  if ($data[1] eq "CDS") {
    $cds_start{$data[0]} = $data[2];
    $cds_end{$data[0]} = $data[3];
  }
  if ($data[1] eq "Transcript") {
    $transcript_start{$data[0]} = $data[2];
    $transcript_end{$data[0]} = $data[3];
    $length{$data[0]} = $transcript_end{$data[0]};
  }

}
close in;

#=head

$row_count = 0;
open (in, "<$treatfile");
while ($line=<in>) {

  $row_count ++;
  if ($row_count%10000 eq 0) {
    print "$row_count\n";
  }

  chomp $line;
  ($rna, $data) = split /\t/, $line;
  @data = split /\;/, $data;
  $sum_treat{$rna} = sum(@data);
}
close in;
print "treat file read in\n";

$row_count = 0;
open (in, "<$controlfile");
while ($line=<in>) {

  $row_count ++;
  if ($row_count%10000 eq 0) {
    print "$row_count\n";
  }

  chomp $line;
  ($rna, $data) = split /\t/, $line;
  @data = split /\;/, $data;
  $sum_control{$rna} = sum(@data);
  if ($length{$rna}) {
    $abt{$rna} = ($sum_treat{$rna}+$sum_control{$rna})/$length{$rna};
#    print "$rna\t$abt{$rna}\n";
    if ($abt{$rna}>$thrd_abundance) {
      $flag{$rna} = 1;
      $count ++;
    }
  }

}
close in;
print "control file read in\n";

print "$count transcripts has abundance >= $thrd_abundance\n";


#=cut


# TREAT and CONTROL score average diagram matrix

$outfile = $filename."_scaled_cds"."$cds_scale"."\_abt$thrd_abundance"."\_treatmatrix.txt";
$outfile2 = $filename."_scaled_cds"."$cds_scale"."\_abt$thrd_abundance"."\_controlmatrix.txt";

#=head

open (out, ">$outfile");
open (out2, ">$outfile2");

for ($j=1; $j<=$total_scale; $j++) {
  push (@js, $j);
}
$js = join "\t", @js;
print out "Base\t$js\n";
print out2 "Base\t$js\n";


open (in, "<$treatfile");
open (in2, "<$controlfile");

%count = ();
%aver = ();


while ($line=<in>) {
  
  $line2 = <in2>;
  $num_row ++;
  if ($num_row%1000 eq 0) {
    print "$num_row\n";
  }
  chomp $line;
  chomp $line2;
  $line =~s/\;/\t/g;
  $line2 =~s/\;/\t/g;
  ($rna, @data) = split /\t/, $line;
  ($rna2, @data2) = split /\t/, $line2;


  for ($i=0; $i<=$#data; $i++) {
    $data[$i] = log($data[$i]+1)/log(2);
    $data2[$i] = log($data2[$i]+1)/log(2);
    $data[$i] = sprintf("%.3f", $data[$i]);
    $data2[$i] = sprintf("%.3f", $data2[$i]);
  }


  $utr5_start = $utr5_start{$rna};
  $utr5_end = $utr5_end{$rna};
  $cds_start = $cds_start{$rna};
  $cds_end = $cds_end{$rna};
  $utr3_start = $utr3_start{$rna};
  $utr3_end = $utr3_end{$rna};


  $utr5_length = $utr5_end{$rna} - $utr5_start{$rna} + 1;
  $cds_length = $cds_end{$rna} - $cds_start{$rna} + 1;
  $utr3_length = $utr3_end{$rna} - $utr3_start{$rna} + 1 - $cut;
  
#    print "$utr5_length; cds $cds_length\t$utr3_length\n";
#    print "$flag{$rna} and $utr5_start{$rna} and $utr3_start{$rna} and ($cds_length >= $cds_scale) and ($utr5_length >= $utr5_scale) and ($utr3_length >= $utr3_scale\n";

  my @data_utr5 = @data[($utr5_start-1)..($utr5_end-1)];
  my @data2_utr5 = @data2[($utr5_start-1)..($utr5_end-1)];
  my @data_cds = @data[($cds_start-1)..($cds_end-1)];
  my @data2_cds = @data2[($cds_start-1)..($cds_end-1)];
  my @data_utr3 = @data[($utr3_start-1)..($utr3_end-1)];
  my @data2_utr3 = @data2[($utr3_start-1)..($utr3_end-1)];

  if ( $flag{$rna} and $utr5_start{$rna} and $utr3_start{$rna} and $cds_start{$rna} ) {
#    print out "$rna";
#    print out2 "$rna";
#    @seqs = split "", $seq{$rna};
     
 #    if ($cds_length >= $cds_scale) {
 
    $utr5_div = $utr5_length/$utr5_scale;
    $cds_div =  $cds_length/$cds_scale;
    $utr3_div = $utr3_length/$utr3_scale;

#    print "$utr5_div, $cds_div, $utr3_div\n";


    print out "$rna";
    print out2 "$rna";

    $window_start = 0;
    for ($j=1; $j<$utr5_scale; $j++) {
      $window_start = ($j-1)*$utr5_div;
      $window_end = $window_start + $utr5_div;
      $int_window_start = int($window_start);
      $int_window_end = int($window_end);
      $sub_length = $int_window_end - $int_window_start + 1;

      @sub_utr5 = @data_utr5[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data_utr5[$int_window_start];
      }
      else {
	$value = sum(@sub_utr5)/$sub_length;
      }
      print out "\t$value";
      $value1 = $value;

 #     print "$j, $window_start, $window_end, $int_window_start, $int_window_end, $data_utr5[$int_window_start], $data_utr5[$int_window_end], $value\n";
   
      @sub_utr5 = @data2_utr5[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data2_utr5[$int_window_start];
      }
      else {
	$value = sum(@sub_utr5)/$sub_length;
      }
      print out2 "\t$value";
      $value2 = $value;

    }
    print out "\t$value1";
    print out2 "\t$value2";

    $window_start = 0;
    for ($j=1; $j<$cds_scale; $j++) {
      $window_start = ($j-1)*$cds_div;
      $window_end = $window_start + $cds_div;
      $int_window_start = int($window_start);
      $int_window_end = int($window_end);
      $sub_length = $int_window_end - $int_window_start + 1;
      
#      print "$j, $window_start, $window_end, $int_window_start, $int_window_end\n";
      
      @sub_cds = @data_cds[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data_cds[$int_window_start];
      }
      else {
	$value = sum(@sub_cds)/$sub_length;
      }
      print out "\t$value";
      $value1 = $value;
      
      @sub_cds = @data2_cds[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data2_cds[$int_window_start];
      }
      else {
	$value = sum(@sub_cds)/$sub_length;
      }
      print out2 "\t$value";
      $value2 = $value;

    }
    
    print out "\t$value1";
    print out2 "\t$value2";
    
    $window_start = 0;
    for ($j=1; $j<=$utr3_scale; $j++) {
      $window_start = ($j-1)*$utr3_div;
      $window_end = $window_start + $utr3_div;
      $int_window_start = int($window_start);
      $int_window_end = int($window_end);
      $sub_length = $int_window_end - $int_window_start + 1;
      
 #     print "$j, $window_start, $window_end, $int_window_start, $int_window_end\n";
      
      @sub_utr3 = @data_utr3[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data_utr3[$int_window_start];
      }
      else {
	$value = sum(@sub_utr3)/$sub_length;
      }
      print out "\t$value";
      
      @sub_utr3 = @data2_utr3[$int_window_start..$int_window_end];
      if ($sub_length eq 1) {
	$value = $data2_utr3[$int_window_start];
      }
      else {
	$value = sum(@sub_utr3)/$sub_length;
      }
      print out2 "\t$value";
      
      
    }
    print out "\n";
    print out2 "\n";

  }
}

close in;
close in2;
close out;
close out2;


#=cut


$averfile = $filename."_scaled_cds"."$cds_scale"."\_abt$thrd_abundance"."\_averageGraph.txt";
open (out, ">$averfile");

%sum = ();
%count = ();
open (in, "<$outfile");
$line = <in>;
print out $line;
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  for ($j=1; $j<=$#data; $j++) {
    if ($noZero) {
      if ($data[$j] ne 0) {
	$sum{$j} = $sum{$j} + $data[$j];
	$count{$j} ++;
      }
    }
    else {
      $sum{$j} = $sum{$j} + $data[$j];
      $count{$j} ++;
    }
  }
}
close in;
print out "treat";
for ($j=1; $j<=$#data; $j++) {
  if ($count{$j}) {
    $aver{$j} = $sum{$j}/$count{$j};
  }
  else {
    $aver{$j} = 0;
  }
  print out "\t$aver{$j}";
}
print out "\n";


%sum = ();
%count = ();
open (in, "<$outfile2");
$line = <in>;
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  for ($j=1; $j<=$#data; $j++) {
    if ($noZero) {
      if ($data[$j] ne 0) {
	$sum{$j} = $sum{$j} + $data[$j];
	$count{$j} ++;
      }
    }
    else {
      $sum{$j} = $sum{$j} + $data[$j];
      $count{$j} ++;
    }
  }
}
close in;
print out "control";
for ($j=1; $j<=$#data; $j++) {
  if ($count{$j}) {
    $aver{$j} = $sum{$j}/$count{$j};
  }
  else {
    $aver{$j} = 0;
  }
  print out "\t$aver{$j}";
}
print out "\n";

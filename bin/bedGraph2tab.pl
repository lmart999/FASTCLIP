#!/usr/bin/perl
#
#######################################################################################
#
# this script is used to find rna per-base value from a .bedGraph file
#
# input files are: annotaion_file, chrsize file, and bedGraph files
# output file is: rna_per-base_signal_tab_file
#
#######################################################################################


print "command line: perl bedGraph2tab.pl annotation_file chrsize_file input.bedGraph output_file\n";

if (!@ARGV) {
  print "command line: perl bedGraph2tab.pl annotation_file chrsize_file input.bedGraph output_file\n";
  exit;
}


$raw_read_length = 1;
$flag_depth_factor = 0;

$location_file = $ARGV[0];
$chrsize = $ARGV[1];
$rnaexp_files = $ARGV[2];
$outfile = $ARGV[3];


# read in chromosome file
open (in, "<$chrsize");
while ($line=<in>) {
  chomp $line;
  ($chr, $size) = split /[\t+\s+]/, $line;
  push (@chrs, $chr);
}
close in;


open (in, "<$location_file");
while ($line=<in>) {
  if (!($line=~/^\#/)) {
    chomp $line;
    @data = split /\t/, $line;
    $transcript_id = $data[1];
    $strand = $data[3];

    $data[9]=~s/\,$//;
    $data[10]=~s/\,$//;

    @exon_starts = split /\,/, $data[9];
    @exon_ends = split /\,/, $data[10];
 
    $effective_length = 0;
    for ($i=0; $i<=$#exon_starts; $i++) {
      $effective_length = $exon_ends[$i] - $exon_starts[$i] + 1 + $effective_length;
    }
    
    $effective_length{$transcript_id} = $effective_length;
    $strand{$transcript_id} = $strand;

    if ( $effective_length eq 0) {
      print "$line\n";
    }
  }
}
close in;


#=head


print "$rnaexp_files\n";  
open (in, "<$rnaexp_files");
  
while ($line=<in>) {
  chomp $line;
  @data = split /[\s+\t+]/, $line;
  $sum = $sum + abs($data[3]*($data[2]-$data[1]));
}

close in;

$actualLen=int($sum/$raw_read_length);
if ($flag_depth_factor) {
  $depth_factor = $actualLen/1000000;
}
else {
  $depth_factor = 1;
}

open (out, ">$outfile");


#@chrs = ("chr1");

foreach $chr (@chrs) {

  undef %flag;
  undef %exp;

  my %flag = ();
  my %exp = ();

  open (in, "<$location_file");
  while ($line=<in>) {
    if (!($line=~/^\#/)) {
      chomp $line;
      @data = split /\t/, $line;
      if ($data[2] eq $chr) {

	$data[9]=~s/\,$//;
	$data[10]=~s/\,$//;

	@exon_starts = split /\,/, $data[9];
	@exon_ends = split /\,/, $data[10];
	for ($i=0; $i<=$#exon_starts; $i++) {
	  for ($j=$exon_starts[$i]; $j<=$exon_ends[$i]; $j++) {
	    $flag{$j} = 1;
	  }
	}
      }
    }
  }
  close in;
  
 
  $file = $rnaexp_files;
  print "$chr\t$file\n";
  
  open (in, "<$file");   
  $count = 0;
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    
    if ($data[0] eq $chr) {
      
      $count ++;
      if ($count%100000 eq 0) {
	print "$chr\t$count\t$line\n";
      }
      
      for ($j=$data[1]+1; $j<=$data[2]; $j++) {
	if ($flag{$j}) {
	  $exp{$j} = abs($data[3])/$depth_factor;
	}
      }
    }
  }
  close in;

  open (in, "<$location_file");
  $count = 0;
  
  while ($line=<in>) {      
    if (!($line=~/^\#/)) {
      chomp $line;
      @data = split /\t/, $line;
    
      $transcript_id = $data[1];
      $strand = $data[3];
      $gene_name = $data[12];
      
      if ($chr eq $data[2]) {
	print out "$transcript_id\t";

	$data[9]=~s/\,$//;
	$data[10]=~s/\,$//;
	
	if ($strand eq "-") {
	  @exon_starts = split /\,/, $data[9];
	  @exon_ends = split /\,/, $data[10];
	  for ($i=$#exon_starts; $i>=0; $i--) {
	    for ($j=$exon_ends[$i]; $j>$exon_starts[$i]; $j--) {
	      if (!$exp{$j}) {
		$exp{$j} = 0;
	      }
	      $value = sprintf("%.3f", $exp{$j});
	      print out "$value\;";
	    }
	  }	  
	  print out "\n";        
	}
	else {	  
	  @exon_starts = split /\,/, $data[9];
	  @exon_ends = split /\,/, $data[10];
	  for ($i=0; $i<=$#exon_starts; $i++) {
	    for ($j=$exon_starts[$i]+1; $j<=$exon_ends[$i]; $j++) {

	      if (!$exp{$j}) {
		$exp{$j} = 0;
	      }
	      $value = sprintf("%.3f", $exp{$j});
	      print out "$value\;";
	    }
	  }
	  print out "\n";
	}
	
      }
    }
  }
  close in;
  
}

close out;

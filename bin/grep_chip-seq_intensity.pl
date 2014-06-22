#!/usr/bin/perl

=head

this script is used to find chip-seq intensity in certain region

input files are: region_location_file in bed format and chip-seq_BedGraph_file
output file is: chip-seq_peak_intensity_file


=cut

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

print "command line: perl grep_chip-seq_intensity.pl location_file file1.wig file2.wig ...\n";


$location_file = $ARGV[0];
shift @ARGV;
@rnaexp_files = @ARGV;


$outfile = "$location_file\_heatmap.txt";
$half_width = 200;
$res = 10;
$raw_read_length = 50;
$low_thrd = 0;


#=head

foreach $file (@rnaexp_files) {
  print "$file\n";

  open (in, "<$file");
  
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s\t]/, $line;
    $sum{$file} = $sum{$file} + $data[3]*($data[2]-$data[1]);
  }
 
  close in;
}


open (out, ">$location_file.total_reads.txt");
foreach $file (@rnaexp_files) {
  my $num_mapped_reads=int($sum{$file}/$raw_read_length);
  print out "$file\t$num_mapped_reads\n";
}
close out;


#=cut


open (in, "<$location_file.total_reads.txt");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  $factor{$data[0]} = $data[1]/1000000;
}
close in;


open (in, "<$location_file");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  $chr = $data[0];
  $flag_chr{$chr} = 1;
  $peak_length = $data[2] - $data[1];
  $num_pdata = int(($peak_length/2 + $half_width)/$res);
}
close in;

open (out, ">$outfile");
print out "peak";
for ($i=0; $i<=$#rnaexp_files; $i++) {
  for ($j=-$num_pdata; $j<$num_pdata; $j++) {
    print out "\t$j";
  }
}
print out "\n";


@chrs = keys %flag_chr;
@chrs = sort {$a cmp $b} @chrs;

#@chrs = ("chr14");

foreach $chr (@chrs) {

  undef %flag;
  undef %exp;

  if ($flag_chr{$chr}) {

    my %flag = ();
    open (in, "<$location_file");
    while ($line=<in>) {
      chomp $line;
      @data = split /\t/, $line;
      if ($data[0] eq $chr) {
	$start = $data[1] - $half_width;
	$end = $data[2] + $half_width;
	for ($i=$start; $i<=$end; $i++) {
	  $flag{$i} = 1;
	}
      }
    }
    close in; 
    
    
    my %exp = ();
    foreach $file (@rnaexp_files) {
      print "$chr\t$file\n";
      
      open (in, "<$file");   
      $count = 0;
      while ($line=<in>) {
	chomp $line;
	@data = split /[\s\t]/, $line;
	
	if ($data[0] eq $chr) {
	  
	  $count ++;
	  if ($count%100000 eq 0) {
	    print "$chr\t$count\t$line\n";
	  }
	  
	  for ($j=$data[1]; $j<$data[2]; $j++) {
	    if ($flag{$j}) {
	      $exp{$j}{$file} = $data[3];
	    }
	  }
	}
      }
      close in;
    }
    
    open (in, "<$location_file");
    $count = 0;
    
    while ($line=<in>) {      
      
      $count ++;
      if ($count%10000 eq 0) {
	print "$count\n";
      }
      
      chomp $line;
      @data = split /\t/, $line;
      $peak_info = "$data[0]:$data[1]-$data[2]";  


      if ($chr eq $data[0]) {
	print out "$peak_info";
	foreach $file (@rnaexp_files) {
	  

	  for ($i=-$num_pdata; $i<$num_pdata; $i++) {
	    $value = 0;

	    for ($j=1; $j<=$res; $j++) {
	      $ith = $data[1] + $res*$i + $j - 1;
	      $value = $value + $exp{$ith}{$file}/$factor{$file};
	    }

	    if ($value >= 1) {
	      $value = log($value)/log(2);
	    }
	    else {
	      $value = 0;
	    }
            
            if ($value<=$low_thrd) {
	      $value = 0;
            }
	
	    print out "\t$value";
	  }
	  
	}	
	print out "\n";
      }
    }
    close in;
  }
}

close out;

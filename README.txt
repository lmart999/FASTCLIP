
---

Steps to start:

(1) Clone the repo
- git clone https://github.com/lmart999/FASTCLIP

(2) Obtain and put a repeat masker file in ~/docs/
- This can be obtained as explained here: http://fantom.gsc.riken.jp/zenbu/wiki/index.php/Uploading_UCSC_repetitive_elements_track
- This is not included in the repo because the file is large.
- The code will, by default, look for a file with name: ~/docs/repeat_masker.bed
- Simply re-name your file to repeat_masker.bed (or modify the code to target your file name).

(3) Obtain and put the hg19 bowtie2 index in ~/docs/hg19/*
- Again, this is not included in the repo because the files are large.
- Files can be obtained here: wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
- Format will be hg19.1.bt2, hg19.3.bt2, etc.

(4) Create a /rawdata directory within the parent (~/rawdata/)
- Move paired raw iCLP .fastq files to /rawdata
- Un-zip and use name convention: <name>_R1.fastq, <name>_R2.fastq.

(5) Create result file directory (~/results/<name>/).
- Output files for <name>_R1/2.fastq will be sent to ~/results/<name>/

---

Dependencies:

(1) Python 2.7 (for CLIPper algorithm)
- https://www.python.org/download/releases/2.7/
(2) iPython 
- http://ipython.org/install.html
(3) iPython notebook (for using the various notebooks provided)
- http://ipython.org/notebook
(4) Matplotlib (plotting)
- http://matplotlib.org/
(5) Pandas (data)
- http://pandas.pydata.org/
(6) Bowtie and bedTools
- http://bowtie-bio.sourceforge.net/index.shtml
- http://bedtools.readthedocs.org/en/latest/
(7) CLIPPER
- https://github.com/YeoLab/clipper/wiki/CLIPper-Home

---

Usage:

(1) Run the main pipeline:
$ ipython fastclip.py.py <name>
(2) The plots will, by defualt, be automatically generated.
(3) The full fastclip pipeline is also provided as a notebook for manual analsis: Fastclip.ipynb.
(4) Several additional notebooks are included see *ipynb.

---

Pipeline steps explained:

(1) Get unzipped reads from the /rawdata directory.
- Format is <NAME>_R< 1 or 2>.fastq.

(2) Trim adapter from the 3' end from the reads.
- RT primer is cleaved, leaving adapters. 
- Remove adapter region from the 3' end of the read.
- The adapter is an input parameter.
- Default is to remove sequnces less than N=33 reads. 
- Q33 specifies the qualtiy score encoding format.

(3) Quality filter.

(4) Remove duplicates.
- This step takes advantage of the fact that 5' end of each read has a random barcode.
- Each initial starting molecule that was RT'd will have a unique barcode.
- Therefore, PCR duplicates are removed by collapsing molecules with identical 5' barcode sequences.

(5) After duplicate removal, remove the 5' barcode sequence.

(6) We then map to a repeat index.
- We use k=1, meaning bt2 will search for 1 distinct, valid alignments for each read.
- This step allows us to both remove reads that are normally blacklisted and also map reads to the repeat index.
- *** The repeat index is derived from < ... >. ***

(7) After mapping, we isolate the 5' position (RT) stop for both positive and negative strand reads.
- This represents the cross-link site in the initial expariment.

(8) For each strand, we merge RT stops between replicates. 
- This means that at RT stop position must be conserved between replicates.
- If conserved, we count the total number of instances of the RT position for both replicates.
- If the total counts exceed a specified threshold, then we record these RT stops.
- Finally, we re-generate a "read" around the RT stop using the passed parameter "expand," 
- A "read" around the RT stop is required for downstream processing.

(9) Reads that do not map to the repeat index are then mapped to hg19.
- The mapped reads are processed by samtools, repeat masker, and blacklist filter.
- As with the repeat index, we then merge RT stops.

(10) Expanded reads from RT stop merging are passed to CLIPper, a peak calling algorithm.
- CLIPper returns a bed-like file format with window coordinates, reads counted per window, etc.
- We use these windows to extract "low FDR" reads from the total set of reads passed to CLIPper.
- The windows provide gene names, which we parse and use to annotate the processed reads.
- We then make bedGraph and BigWig files from this complete pool of "low FDR" reads, allowing easy visualization. 

(11) Partition "low FDR" reads by gene type.
- We partition the gene names recoved from CLIPper using ENSEMBL annotation of gene name by RNA type.
- Once this is done, we also split the "low FDR" reads recovered from CLIPper by type using the gene name.
- Protein coding and lincRNA genes can be embedded snoRNAs or miRNAs that make this more challenging.
- In turn, we re-generate the initial RT stop and intersect this with two different filters.
- One filter is a snoRNA mask and the other filter is a miR mask. 
- *** Both are derived from < ... >. ***
- These masks allow us to remove all "protein coding" RT stops that fall within annotated sno/mi-RNA regions.

(12) Quantification of reads per gene.
- For each gene type, we quantify the number of reads per gene.
- For all but snoRNAs, this is computed using the bed files obtained above.
- For snoRNAs, we intersect the initial pool of "low FDR" reads with custom annotation file. 
- Collectivly, this gives us reads per gene for each gene type.
- All are based upon ENSEMBL annotation except for the snoRNAs.

(13) Summary of RT stop intensity around CLIPper cluster centers.
- We generate a bed file of cluster center positions using the CLIPper cluster file output.
- We use a custom perl script that generates a heatmap of RT stop intensty per cluster.
- This allows us to later visualize the distribution of RT stops per cluster.

(14) Parition protein coding reads by UTR.
- We intersect sno/mi-RNA filtered reads with ENSEMBL-derived UTR coordinates. 
- We perform this such that each read assignment is mutually exclusive.
- This only isolates reads that fall within each UTR type.
- Similarly, we use a custom perl script to generate a matrix of read intensity per gene.
- This provides a complete binding profile per gene.

(15) Partition reads by ncRNA binding region.
- For non-coding RNAs, we simply annotate reads with the start and stop position for each ncRNA.
- This allows us to determine the position of each RT stop with respect to the full length of the gene.

(16) Partition repeat-mapped RT stops by region.
- The repeat RNA mapped RT stops are paritioned using the repeat custom index annotation. 
- As with the ncRNAs, this is later used for visualization.

(16) Figure 1 visualizes the some of the relevant summary data.
- It includes a read count summary per pipeline step. The source data is: PlotData_ReadsPerPipeFile
- It includes a pie chart of UTR binding.
- This uses reads obtrained from intersection with ENSEMBL-derived UTR coordinates. 
- The source data is: PlotData_ReadsPerGene_*UTR or CDS
- It also includes a bar graph of gene count per RNA type.
- The source data is: PlotData_ReadAndGeneCountsPerGenetype
 
(17) Figure 2 provides a richer summary of the UTR data.
- The upper panel is an aggregate trace of binding derived from a custom perl script.
- The lower panels provide a heatmap of binding intensity for gene exclusivly bound in each UTR or CDS.
- This allows us to isolate genes with exclusive UTR,CDS,or intronic binding.
- The source data is: PlotData_ExclusiveBound_*

(18) Figure 3 and 4 provides coverage histograms of binding across each repeat RNA and rRNA, respectivly. 
- Source data: PlotData_RepeatRNAHist_*

(19) Figure 5 provides a summary of snoRNA binding data.
- The pie chart provides a summary of reads per snoRNA type.
- These are complimented by histograms of RT stop position within the snoRNA gene body.

(20) Figure 6 provides histograms of RT stop position within gene body for all remaining ncRNA types.

---

Debugging:

(1) Mapping
- Ensure that bowtie2 is in the $PATH and executable.

(2) CLIPPER
- Uses Python 2.7. 

(3) Any scrip in /bin 
- The provided BedGraphToBigWig is built for Linux OS.
- This, and related scripts, may be downloaded for other platforms: 
http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.i386/

(4) If there is a problem parsing clusters from CLIPper, consider the version used and update the CLIPPERoutNameDelim='_' parameter.
- The older version of CLIPper results <name>_<clusterNum>_<readPerCluster>
- The newer version has name.val__<clusterNum>_<readPerCluster>
- Therefore, for newer version set CLIPPERoutNameDelim='.'

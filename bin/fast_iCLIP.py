

#!/usr/bin/env python
import os
import cmath
import math
import sys
import numpy as np
import glob 
import subprocess
import re
# from matplotlib_venn import venn2
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
from optparse import OptionParser

def getFastq(infilepath,sampleName):
	# Useage: Get the path to a pair of zipped files (either replicates or paired-end reads).
	# Input: Path to raw data files ans sample name.
	# Output: The path to each read.
	read1=infilepath+sampleName+'_R1.fastq'
	read2=infilepath+sampleName+'_R2.fastq'
	return(read1,read2)
 
def moveFiles(reads,outfilepath):
	# Useage: Copy a list of files of a given destination and updates their path.
	# Input: List of files and a new destination directory.
	# Output: Updated list of file paths.
	try:
		for read in reads:
			proc = subprocess.Popen(['cp',read,outfilepath], shell=False, stderr=subprocess.PIPE)
			proc.communicate()
		return changePath(reads,outfilepath) 
	except:
		print "Error with moving files."

def changePath(inpaths,outpath):
	# Useage: Change the path of a file.
	# Input: List of path for files and new path for all.
	# Output: List of new paths.
	pathNames=[]
	try:
		for inpath in inpaths:
			head, tail = os.path.split(inpath)
			newfilepath=outpath+tail
			pathNames = pathNames+[newfilepath]
		return pathNames
	except:
		print "Problem changing file handle."

def changeHandle(reads,newHandle):
	# Useage: Change the handle of a file.
	# Input: File of type <name>.<handle> 
	# Output: List of <name>.<new handle>
	readNames=[]
	try:
		for read in reads:
			head, tail = os.path.split(read)
			newfilepath=head+'/'+tail.split('.')[0]+newHandle
			readNames = readNames+[newfilepath]
		return readNames
	except:
		print "Problem changing file handle."

def changePath(inpaths,outpath):
	# Useage: Change the path of a file.
	# Input: List of path for files and new path for all.
	# Output: List of new paths.
	pathNames=[]
	try:
		for inpath in inpaths:
			head, tail = os.path.split(inpath)
			newfilepath=outpath+tail
			pathNames = pathNames+[newfilepath]
		return pathNames
	except:
		print "Problem changing file handle."
	
def countFiles(unzippedreads):
	# Useage: Counts the number of lines in a specified file.
	# Input: List of fastq files.
	# Output: Prints the line count to standard out.
	try:
		for read in unzippedreads:
			print read
			process=subprocess.Popen(['wc','-l',read], shell=False, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			print process.communicate()
	except:
		print "Problem with counting file."

def modifyName(filepath,newTag):
	# Useage: Modifies the filepath name. 
	# Input: File path of format <path>/<name>.fastq and a string to add to the name.
	# Output: Returns the modified path of type <old path>_<new modifier>.fastq
	try:
		head, tail = os.path.split(filepath)
		oldname = tail.split('.')[0]
		newName = head+"/"+oldname+"_"+newTag
		return newName
	except:
		print "Problem with modifying file name."

def trimReads5p(unzippedreads,n):
	# Useage: Trims a specified number of bases from the 5' end of each read.
	# Input: List of fastq files.
	# Output: List of 5p trimmed files.
	trimparam='-f'+str(n)
	trimmedReads=[]
	print "Performing raw data filtering.."
	logOpen.write("Performing 5p trimming...\n")
	try:
		for inread in unzippedreads:
			outread = modifyName(inread,"5ptrimmed.fastq")
			# -Q33 indicates Illumina quality score encoding
			process=subprocess.Popen(['fastx_trimmer', trimparam, '-Q33', '-i', inread,'-o',outread], shell=False,stderr=subprocess.PIPE)
			# Communicate the process so that the function waits to finish before exiting.
			process.communicate()
			trimmedReads=trimmedReads+[outread]
		return trimmedReads
	except:
		logOpen.write("Problem with 5p trimming.\n")
		print "Problem with 5p trimming."

def trimReads3p(unzippedreads,adapter3p):
	# Useage: Trims a specified adapter sequence from the 3p end of the reads.
	# Input: List of (5' trimmed) fastq files.
	# Output: List of 3p trimmed files.
	trimparam='-a'+adapter3p
	trimmedReads=[]
	logOpen.write("Performing 3p trimming...\n")
	try:
		for inread in unzippedreads:
			outread = modifyName(inread,"3ptrimmed.fastq")
			# Parameters:
			# -n: keep sequences with unknown (N) nucleotides
			# -D: DEBUG output
			# -l: Discard sequences shorter than N nucleotides
			# -i $FASTQ.fastq -o
			process=subprocess.Popen(['fastx_clipper', trimparam, '-n', '-l33', '-Q33', '-i', inread,'-o',outread], shell=False,stderr=subprocess.PIPE)
			process.communicate()
			trimmedReads=trimmedReads+[outread]
		return trimmedReads
	except:
		logOpen.write("Problem with 3p trimming.\n")
		print "Problem with 3p trimming."

def qualityFilter(unzippedreads,q,p):
	# Useage: Filters reads based upon quality score.
	# Input: List of fastq file names as well as the quality paramters p and q.
	# Output: List of modified fastq file names.
	qualityparam='-q'+str(q)
	percentrageparam='-p'+str(p)
	filteredReads=[]
	logOpen.write("Performing quality filtering...\n")
	try:
		for inread in unzippedreads:
			outread = modifyName(inread,"filter.fastq")
			# Parameters:
   			# q: Minimum quality score to keep.
	   		# p: Minimum percent of bases that must have [-q] quality.
			process=subprocess.Popen(['fastq_quality_filter', qualityparam, percentrageparam,'-Q33', '-i', inread,'-o',outread], shell=False,stderr=subprocess.PIPE)
			process.communicate()
			filteredReads=filteredReads+[outread]
		return filteredReads
	except:
		logOpen.write("Problem with quality filter.\n")
		print "Problem with quality filter."

def dupRemoval(unzippedreads):
	# Useage: Removes duplicate reads.
	# Input: List of fastq file names.
	# Output: List of reads in FASTA format.
	filteredReads=[]
	logOpen.write("Performing duplicate removal...\n")
	print "Performing duplicate removal..."
	try:
		for inread in unzippedreads:
			outread = modifyName(inread,"nodupe.fasta")
			process=subprocess.Popen(['fastx_collapser','-Q33', '-i', inread,'-o',outread], shell=False,stderr=subprocess.PIPE)
			process.communicate()
			filteredReads=filteredReads+[outread]
		return filteredReads
	except:
		logOpen.write("Problem with duplicate removal.\n")
		print "Problem with duplicate removal."

def fastaTofastq(fastaIN):
	# Usage: Convert fasta to fastq
	# Input: List of fasta files
	# Output: List of fastq files
	program =  os.getcwd() + '/bin/fasta_to_fastq.pl'
	fastqFiles=[]
	print "Performing fasta conversion..."
	try:
		for fasta in fastaIN:
			fastqOut = fasta.replace('.fasta', '.fastq')
			outfh = open(fastqOut, 'w')
			proc = subprocess.Popen(['perl',program,fasta],stdout=outfh)
			proc.communicate()
			fastqFiles=fastqFiles+[fastqOut]
		return fastqFiles
	except:
		print "Error with fastq file conversion."

def runBowtie(fastqFiles):
	# Useage: Short read mapping to reference (hg19).
	# Input: Fastq files of replicates (not paired end).
	# Output: Path to samfile for each read.
	program = 'bowtie2'
	mappedReads=[]
	unMappedReads=[]
	print "Performing Bowtie..."
	logOpen.write("Performing Bowtie...\n")
	# Parameters
	# -m : Suppress all alignments for a particular read or pair if more than <int> reportable alignments exist for it: -m 1
	# -v : Alignments may have no more than `V` mismatches: -v2
	try:
		for infastq in fastqFiles:
			print "Input file:"
			print infastq 
			# Process the genome index
			if genomeDict[genomeIndex] == 1:
				index = os.getcwd() + '/docs/hg19/hg19'
				outfile = modifyName(infastq,"mapped.sam")
				unmapped = modifyName(infastq,"notMappedToHg19.fastq")
			elif genomeDict[genomeIndex] == 2:
				index = os.getcwd() + '/docs/jfh1/jfh1'
			elif genomeDict[genomeIndex] == 3:
				index = os.getcwd() + '/docs/h77/h77'
			elif genomeDict[genomeIndex] == 4:
				index = os.getcwd() + '/docs/mm9/mm9'
			elif genomeDict[genomeIndex] == 5:
				index = os.getcwd() + '/docs/repeat/rep'
				outfile = modifyName(infastq,"mappedToRepeatRNA.sam")
				unmapped = modifyName(infastq,"notMappedToRepeat.fastq")
			print 'Genome index:'
			print index
			print "Output file (mapped):"
			print outfile
			print "Output file (unmapped)"
			print unmapped
			proc = subprocess.Popen([program,'-x', index,'-U',infastq,'--un',unmapped,'-S',outfile,'2>>%s'%logFile],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			proc.communicate()
			mappedReads = mappedReads + [outfile]
			unMappedReads = unMappedReads + [unmapped]
		return (mappedReads,unMappedReads)
	except:
		logOpen.write("Error with Bowtie.\n")
		print "Error with Bowtie"

def runSamtools(samfiles):
	# Useage: Samfile processing.
	# Input: Sam files from Bowtie mapping.
	# Output: Duplicate removed, sorted bedFiles.
	program = 'samtools'
	program2 = 'bamToBed'
	outBedFiles=[]
	logOpen.write("Performing Samtools...\n")
	try:
		
		for samfile in samfiles:

			if genomeIndex == 'rep':

				# Convert to bamfile
				bamfile = samfile.replace('.sam', '.bam')  
				proc = subprocess.Popen( [program, 'view', '-bS', '-o', bamfile, samfile])
				proc.communicate()

				# Sort the bamfile and note that samtools sort adds the .bam handle
				bamfile_sort = bamfile.replace('.bam', '_sorted') 
				proc2 = subprocess.Popen( [program, 'sort', bamfile, bamfile_sort])
				proc2.communicate()
				
				# Convert to bedFile
				bedFile = bamfile_sort.replace('_sorted', '.bed')
				outfh = open(bedFile, 'w')
				proc3 = subprocess.Popen( [program2,'-i', bamfile_sort+'.bam'],stdout=outfh)
				proc3.communicate()

			else:

				# Convert to bamfile
				bamfile = samfile.replace('.sam', '.bam')  
				proc = subprocess.Popen( [program, 'view', '-bS', '-o', bamfile, samfile])
				proc.communicate()

				# Sort the bamfile and note that samtools sort adds the .bam handle
				bamfile_sort = bamfile.replace('.bam', '_sorted') 
				proc2 = subprocess.Popen( [program, 'sort', bamfile, bamfile_sort])
				proc2.communicate()

				# Remove duplicates
				bamfile_nodupes = bamfile_sort.replace('_sorted', '_nodupes.bam') 
				proc3 = subprocess.Popen( [program, 'rmdup','-s', bamfile_sort+'.bam', bamfile_nodupes])
				proc3.communicate()

				# Convert to bedFile
				bedFile = bamfile_nodupes.replace('_nodupes.bam', '_nodupes.bed')
				outfh = open(bedFile, 'w')
				proc4 = subprocess.Popen( [program2,'-i', bamfile_nodupes],stdout=outfh)
				proc4.communicate()

			outBedFiles=outBedFiles+[bedFile]

		return outBedFiles

	except:
		logOpen.write("Error with Samtools.\n")
		print "Error with Samtools"

def seperateStrands(mappedReads):
	# Useage: Seperate positive and negative strands.
	# Input: Paths to two bed files from Samtools.
	# Output: Paths to bed files isolated by strand.
	logOpen.write("Performing isolation of reads by strand...\n")
	try:
		# Create list for storing file names
		negativeStrand=[]
		positiveStrand=[]
		# For each file in the input list
		for mapFile in mappedReads:
			# Open the file
			with open(mapFile, 'r') as infile:
				# Create new file handles
				neg_strand=modifyName(mapFile,'neg.bed')
				pos_strand=modifyName(mapFile,'pos.bed')	
				# Open new files
				neg = open(neg_strand, 'w')
				pos = open(pos_strand, 'w')
				negativeStrand=negativeStrand+[neg_strand]
				positiveStrand=positiveStrand+[pos_strand]
				# Read one line at the time to memory, and write to outputs.
				for line in infile:	
					# Sort read based upon strand, which is field six (index=5) 
					if str(line.strip().split('\t')[5]) == '-':
						neg.write(line)
					elif str(line.strip().split('\t')[5]) == '+':
						pos.write(line)
		return (negativeStrand,positiveStrand)

	except:	
		logOpen.write("Error with seperating strands.\n")
		print "Error with seperating strands."	

def modifyNegativeStrand(negativeStrandReads):
	# Useage: For negative stranded reads, ensure 5' position (RT stop) is listed first.
	# Input: Bed file paths to all negative stranded.
	# Output: Paths to modified bed files.
	logOpen.write("Modifying the negative strand reads.\n")
	negativeStrandEdit=[]
	try:
		# For each file in the input list
		for negativeRead in negativeStrandReads:
			# Outpit file name
			neg_strand_edited = modifyName(negativeRead,'edit.bed')
			negativeStrandEdit=negativeStrandEdit+[neg_strand_edited]
			# Open new files
			neg_edit = open(neg_strand_edited, 'w')
			with open(negativeRead, 'r') as infile:
				for line in infile:	
					chrom,start,end,name,quality,strand=line.strip().split('\t')
					# For negative stranded reads, invert so that 5' position is listed first (3' position is an arbitrary 30 bases beyond)
					neg_edit.write('\t' .join((chrom,end,str(int(end)+30),name,quality,strand,'\n')))

		return negativeStrandEdit

	except:
		logOpen.write("Error with correcting the negative strand.\n")
		print "Error with correcting the negative strand."	

def isolate5prime(strandedReads):
	# Useage: Isolate only the Chr, 5' position (RT stop), and strand.
	# Input: Bed file paths to strand seperated reads.
	# Output: Paths to 5' isolated reads.
	logOpen.write("Isolating RT stops.\n")
	try:
		# For each file in the input list
		RTstops=[]
		for reads in strandedReads:
			# Outpit file name
			RTstop = modifyName(reads,'RTstop.bed')
			# Open new files
			f = open(RTstop, 'w')
			with open(reads, 'r') as infile:
				RTstops=RTstops+[RTstop]
				for line in infile:	
					chrom,start,end,name,quality,strand=line.strip().split('\t')
					
					f.write('\t' .join((chrom,start,strand,'\n')))

		return RTstops

	except:
		logOpen.write("Error with isolating RT stops.\n")
		print "Error with isolating RT stop."

def mergeRT(RTstops,outfilename):
	# Useage: Merge RT stops between replicates.
	# Input: Paths to RT stop files (stranded reads) and output filename.
	# Output: Nothing (outfile name is specified in input)
	logOpen.write("Merging RT stops.\n")
	try:
		# Create object for storing dictionaries
		store=[0,0]
		# Flag for storing each dictionary
		i=0
		# Iterare through each file in input
		for RT in RTstops:
			# Create a dictionary with each value initialized to zero
			d = defaultdict(int)
			# Open the file 
			with open(RT, 'r') as infile:		
				for line in infile:	
					if line in d:
						d[line] += 1
					else:
						# Make each RT stop a unique key in the dictionary
						d[line] += 1
			# Store the dictionary
			store[i]=d
			i += 1

		rt_rep1=[k for k in store[0].iteritems()]
		rt_rep2=[k for k in store[1].iteritems()]

		# Open output files
		f = open(outfilename, 'w')

		# Iderate through each RT stop in dictionary 1	
		for key in store[0]:
			# Check if same RT stop is in dictionary 2
			name=store[1].get(key,None)
			# If so, then the read is preserved in both files
			if name:
				chrom,start,strand=key.strip().split('\t')
				# Make sure the start of the read is greater than 15 bases from end of the chrom
				if int(start)>15:
					# Create a new read centered +/- 15 bases around the RT stop coordinate
					read='\t' .join((chrom,str(int(start)-15),str(int(start)+15),'CLIPread','255',strand,'\n'))
				else:
					read='\t' .join((chrom,str(int(start)),str(int(start)+15),'CLIPread','255',strand,'\n'))
				# Write the read to the output, and repeat for the total number of instances the RT stop appears
				f.write(read*(store[0][key]+store[1][key]))

	except:

		logOpen.write("Error with merging RT stops\n")
		print "Error with merging RT stops."

def fileCat(destinationFile,fileList):
	# Useage: Concatenate two files.
	# Input: Output file path, as well as a list input files.
	# Output: Nothing (outfile name is specified in input).
	logOpen.write("Concatening files.\n")
	try:
		f = open(destinationFile, "w")
		for tempfile in fileList:
			# Read each file into the destrination file
			readfile = open(tempfile, "r")
			f.write(readfile.read())
			readfile.close()
		f.close()

	except:

		logOpen.write("Error with file concatenation.\n")
		print "Error with file concatenation."

def runCLIPPER(RTclusterfile):
	# Useage: Process the mergedRT file and pass through CLIPper FDR script.
	# Input: Merged RT file.
	# Output: CLIPper input (.bed) file and output file.
	program = 'bedToBam'
	genomeFile = os.getcwd()+'/docs/human.hg19.genome'
	program2 = 'samtools'
	program3 = 'bamToBed'
	program4 = 'clipper'
	print "Running CLIPper..."
	logOpen.write("Running CLIPper...\n")

	try:
		# Create and open bamfile
		bamfile = RTclusterfile.replace('.bed', '.bam')  
		outfh = open(bamfile, 'w')
		proc = subprocess.Popen([program, '-i', RTclusterfile,'-g',genomeFile],stdout=outfh)
		proc.communicate()

		bamfile_sort = bamfile.replace('.bam', '.srt')
		proc2 = subprocess.Popen([program2, 'sort', bamfile, bamfile_sort])
		proc2.communicate()

		bamfile_sorted=bamfile_sort+'.bam'
		mapStats = bamfile_sorted.replace('.srt.bam', '.mapStats.txt') 
		outfh = open(mapStats, 'w')
		proc3 = subprocess.Popen([program2, 'flagstat', bamfile_sorted],stdout=outfh)
		proc3.communicate()

		proc4 = subprocess.Popen([program2, 'index', bamfile_sorted])
		proc4.communicate()

		CLIPPERin = bamfile_sorted.replace('.srt.bam', '_CLIPPERin.bed') 
		outfh = open(CLIPPERin, 'w')
		proc5 = subprocess.Popen([program3, '-i', bamfile_sorted],stdout=outfh)
		proc5.communicate()

		CLIPPERout = CLIPPERin.replace('_CLIPPERin.bed', '_CLIP_clusters') 
		proc6 = subprocess.Popen([program4, '--bam', bamfile_sorted,'-shg19','--outfile=%s'%CLIPPERout],)
		proc6.communicate()
		outfh.close()

		return (CLIPPERin,CLIPPERout)

	except:
		logOpen.write("Error with running CLIPper.\n")
		print "Error with running CLIPper."

def modCLIPPERout(CLIPPERin,CLIPPERout):
	# Usage: Process the CLIPper output and isolate lowFDR reads based upon CLIPper windows.
	# Input: .bed file passed into CLIPper and the CLIPper windows file.
	# Output: Low FDR reads recovered using the CLIPer windows file, genes per cluster, gene list of CLIPper clusters, and CLIPper windows as .bed.
	program = 'intersectBed'
	logOpen.write("Processing CLIPper output...\n")
	try:
		# Output format from CLIPper will be: <same_name>_CLIP_clusters
		CLIPperOutBed=CLIPPERout+'.bed'
		CLIPpeReadsPerCluster=CLIPPERout+'.readsPerCluster'
		CLIPpeGeneList=CLIPPERout+'.geneNames'
		
		# Open new files
		f = open(CLIPperOutBed, 'w')
		g = open(CLIPpeReadsPerCluster, 'w')
		h = open(CLIPpeGeneList, 'w')
		with open(CLIPPERout, 'r') as infile:
			# For each CLIPper window.
			for line in infile:	
				try:
					# Version of CLIPper used here includes a header that cannot be parsed. Handle this.
					chrom,start,end,name,stats,strand,start_2,end_2 = line.strip().split('\t')
					# Ensembl genes are parsed with <name>.<value>_<cluster>_<count>
					readPerCluster=(name.strip().split('.')[1]).split('_')[2]
					geneName=(name.strip().split('.')[0])
					# Re-write the CLIPper windows file
					f.write('\t' .join((chrom,start,end,name,stats,strand,'\n')))
					g.write((readPerCluster+'\n'))
					h.write((geneName+'\n'))
				except:
					logOpen.write("Problem with CLIPper ID. Continue reading windows file...\n")
		f.close()
		g.close()
		h.close()

		# File name for low FDR reads
		CLIPPERlowFDR = CLIPperOutBed.replace('.bed', '_lowFDRreads.bed')
		# Intersect input reads with the CLIPper windows
		outfh = open(CLIPPERlowFDR, 'w')
		# Note that -u and -wb are mutually exclusive, but enforce strandedness
		proc = subprocess.Popen([program,'-a', CLIPPERin, '-b', CLIPperOutBed,'-wb','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()

		return (CLIPPERlowFDR,CLIPpeReadsPerCluster,CLIPpeGeneList,CLIPperOutBed)

	except:

		logOpen.write("Problem obtaining lowFDR reads.\n")
		print "Problem obtaining lowFDR reads."

def compareLists(list1,list2,outname):
	# Usage: Compare gene lists and output matches to the file. 
	# Input: Two gene lists.
	# Output: Path file containing the matching genes.
	logOpen.write("Comparing gene lists...\n")
	try:
		# Set comparison, resulting in shared genes 
		f = open(list1, 'r')
		g = open(list2, 'r')
	 	content1 = set(f.readlines())
	 	content2 = set(g.readlines())
	 	commonGenes = content1 & content2

	 	# Write shared genes to an output file
	 	geneCategory=outname.split('.')[1]
	 	outputName=outfilepath+'clipGenes_'+geneCategory
	 	outfh = open(outputName, 'w')
	 	for gene in commonGenes:
	 		outfh.write(gene)
	 	outfh.close()

		return outputName

	except:

		logOpen.write("Problem comparing two sets of genes.\n")
		print "Problem comparing two sets of genes."

def getLowFDRGeneTypes(CLIPpeGeneList):
	# Usage: Get all genes listed under each type, compare to CLIPer targets.
	# Input: .bed file passed into CLIPper and the CLIPper windows file.
	# Output: Path to file containing all CLIPper genes of each type.
	logOpen.write("Grabbing genes of each type...\n")
	# try:
	readListByGeneType=[]
	# Iterate through all gene types
	for geneType in os.listdir(os.getcwd() + '/docs/genes_types'):
		# Paths to gene lists and output directory
		genepath=os.getcwd() + '/docs/genes_types/' + geneType
		# Compare list of lowFDR CLIP gene names and gene names by gene type
		lowFDRreadlist=compareLists(CLIPpeGeneList,genepath,geneType)
		# Update the list of paths to the resulting files
		readListByGeneType=readListByGeneType+[lowFDRreadlist]
	return readListByGeneType

def grep(pattern,filein):
	# Usage: Open a file and search all lines for a pattern.
	# Input: Pattern to search (gene name) and file name.
	# Output: List of lines in a file that have the pattern.
	r = []
	filein_open = open(filein, 'r')
	for line in filein_open:
	 	if re.search(pattern,line):
	 		r.append(line)
	filein_open.close()
	return r

def getLowFDRReadTypes(CLIPPERlowFDR,pathToGeneLists):
	# Usage: Given a list of genes, return all reads for the associated genes.
	# Input: Gene list and the path to lowFDR read file.
	# Output: List of reads assocaited with the given genes.
	print "Grabbing low FDR read by gene type..."
	logOpen.write("Grabbing low FDR reads... \n")
	lowFDRgenelist=[]
	try:
		for path in pathToGeneLists:
			# File path to which low FDR reads of each type will be written
			outfile=path+'_LowFDRreads.bed'
			outfh = open(outfile, 'w')
			# searchfile = open(CLIPPERlowFDR, 'r')
			# Open the list of genes for each type
			with open(path, 'r') as infile:
				# For each gene, loopup the corresponding reads
				for geneName in infile:	
					# Grep it to the output 
					store=grep(geneName.strip(), CLIPPERlowFDR)
					# If NOT empty, then write to output
					if store:				
						outfh.write(''.join(store))
			lowFDRgenelist=lowFDRgenelist+[outfile]
			outfh.close()
		return lowFDRgenelist
	except:
		logOpen.write("Problem isolating low FDR reads by type.\n")
		print "Problem isolating low FDR reads by type."					

def lineCountForLog(infileList,outFileName):
	# Usage: Obtain line count for all files in alist.
	# Input: File list and output file path.
	# Output: Path to count for all files in the input list.
	try:
		outfh = open(outFileName, 'w')
		for infile in infileList:
			with open(infile) as fin:
				lines = sum(1 for line in fin)
			outfh.write(str(infile)+'\t'+str(lines)+'\n')
	except:
		"Print error with line count"

def lineCount(filename):
	# Usage: Get and return line count for a file.
	# Input: File 
	# Outline: Lines
	print "Counting lines for ..."
	print filename
	i=0
	with open(filename) as f:
		for i,l in enumerate(f):
			pass
	print "Lines counted:"
	print i
	return i+1

def cleanBedFile(inBed):
	# Usage: Sort and recover only first 6 fields from a bed file.
	# Input: BedFile.
	# Output: Sorted bedFile with correct number of fields.
	program='sortBed'
	try:
		# Make sure bedfile only has 5 fields
		CLIPperOutBed=inBed.replace('.bed','_cleaned.bed')	
		sortedBed=CLIPperOutBed.replace('_cleaned.bed','_cleaned_sorted.bed')
		
		# Open new files
		f = open(CLIPperOutBed, 'w')
		with open(inBed, 'r') as infile:
			for line in infile:	
				elementList = line.strip().split('\t')
				# Re-write the CLIPper windows file
				f.write('\t' .join((elementList[0],elementList[1],elementList[2],elementList[3],elementList[4],elementList[5],'\n')))
		f.close()

		# Sort the resulting bedFile
		outfh = open(sortedBed, 'w')
		proc = subprocess.Popen([program, '-i', CLIPperOutBed],stdout=outfh)
		proc.communicate()
		outfh.close()

		return sortedBed

	except:
		print "Error cleaning file."

def makeBedGraph(lowFDRreads):
	# Usage: From a bedFile, generate a bedGraph and bigWig.
	# Input: BedFile.
	# Output: BedGraph file.
	program = 'genomeCoverageBed'
	program2 = os.getcwd() + '/bin/bedGraphToBigWig'
	# Note include this reference to account for the fact that non-traditional chroms (e.g., chrUn_gl000220) can appear in the data.
	# Without this, bedGraph creation will throw an error.
	sizesFile = os.getcwd()+'/docs/human.hg19.genome'
	sizesFile2 = os.getcwd()+'/docs/hg19.sizes'

	try:
		cleanBed = cleanBedFile(lowFDRreads)
		outname = cleanBed.replace('.bed','.bedgraph')
		outname2 = cleanBed.replace('.bed','.bw')

		outfh = open(outname, 'w')
		proc = subprocess.Popen([program, '-bg','-split','-i',cleanBed,'-g',sizesFile],stdout=outfh)
		proc.communicate()
 
		outfh2 = open(outname2, 'w')
		proc2 = subprocess.Popen([program2,outname,sizesFile,outname2],stdout=subprocess.PIPE)
		proc2.communicate()
		return outname

	except:
		print "Problem making bedGraph."

def makeClusterCenter(windowsFile):
	# Usage: Generate a file of cluster centers.
	# Input: Raw CLIPper output file.
	# Output: File with coordinates for the center of each CLIPper cluster.
	try:
		cleanBed = cleanBedFile(windowsFile)
		centers=cleanBed.replace('.bed','.clusterCenter')
		# Open new files
		f = open(centers, 'w')
		with open(cleanBed, 'r') as infile:
			for line in infile:
				# Each line from the CLIPper cluster file
				elementList = line.strip().split('\t')
				# Lower coordinate in order to correct for differences in strand
				minCoordinate=min(int(elementList[1]), int(elementList[2]))
				diff=abs(int((int(elementList[1])-int(elementList[2]))/2))
				# Write the center of each window
				f.write(elementList[0] + '\t' + str(minCoordinate+diff) + '\t' + str(minCoordinate+diff+1) + '\n')
		f.close()

		return centers

	except:
		print "Problem making the cluster center file."

def getClusterIntensity(bedGraph,centerCoordinates):
	# Usage: Generate a matrix of read itensity values around CLIPper cluster center.
	# Input: BedGraph and cluster center file.
	# Output: Generates a matrix, which is passed into R.
	program = os.getcwd() + '/bin/grep_chip-seq_intensity.pl'
	program2 = 'wait'
	logOpen.write("Generating cluster intensity... \n")
	try:
		proc = subprocess.Popen(['perl',program, centerCoordinates, bedGraph],)
		proc.communicate()

		print "Waiting for Cluster Intensity file completion..."
		proc2 = subprocess.Popen(program2,shell=True)
		proc2.communicate()

	except:
		logOpen.write("Problem with generating cluster intensity.\n")
		print "Problem with get reads around cluster centers."

def makeTab(bedGraph):
	program = os.getcwd() + '/bin/bedGraph2tab.pl'
	program2 = 'wait'
	genesFile = os.getcwd() + '/docs/hg19_ensembl_genes.txt'
	sizesFile = os.getcwd() + '/docs/hg19.sizes'
	try:
		outfile=bedGraph.replace('.bedgraph','.tab')
		print "Waiting for Tabfile completion..."
		proc = subprocess.Popen(['perl',program,genesFile,sizesFile,bedGraph,outfile],)
		proc.communicate()

		proc2 = subprocess.Popen(program2,shell=True)
		proc2.communicate()

		return outfile

	except:
		logOpen.write("Problem making tab.\n")
		print "Problem with making the tab file."

def makeAvgGraph(bedGraph):
	# Usage: Generate a matrix of read itensity values across gene body.
	# Input: BedGraph.
	# Output: Generates two matricies, which are passed into R.
	program= os.getcwd() + '/bin/averageGraph_scaled_tab.pl'
	program2 = 'wait'
	utrFile = os.getcwd() + '/docs/hg19_ensembl_UTR_annotation.txt'
	print "Make average graph..."
	logOpen.write("Generating average graph... \n")
	try:
		tabFile=makeTab(bedGraph)
		outhandle=tabFile.replace('.tab','_UTRs')
		proc = subprocess.Popen(['perl',program,utrFile,tabFile,tabFile,outhandle],)
		proc.communicate()

		# Perl send this to a background process, so wait for completion before going ahead.
		print "Waiting for AverageGraph completion..."
		proc2 = subprocess.Popen(program2,shell=True)
		proc2.communicate()

	except:
		logOpen.write("Problem with making average graph.\n")
		print "Problem with making average graph."

def extractClusters(geneList,allClusters):
	# Usage: Extract a set of CLIPper clusters based on gene name
	# Input: Gene list
	# Output: Clusters with those genes
	try:
		# Create the outfile 
		extractedClusters = geneList+'clusters'
		outfh = open(extractedClusters, 'w')
		# Iterate through each gene name
		namesToQuery = np.genfromtxt(geneList,usecols=(0,),delimiter='\t',dtype='string')
		for name in namesToQuery:
			# Grep it to the output 
			store=grep(name.strip(),allClusters)
			# If NOT empty, then write to output
			if store:				
				outfh.write(''.join(store))
		outfh.close()
		return extractedClusters
	except:
		print "Problem extracting clusters for gene list."

def runBlacklistRegions(mappedReads):
	# Usage: Remove blacklisted regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools).
	# Output: Bedfile with blacklisted regions removed.
	program = 'intersectBed'
	blacklistregions = os.getcwd() + '/docs/wgEncodeDukeMapabilityRegionsExcludable.bed'
	blackListed=[]
	try:
		for bedIn in mappedReads:
			# File name for low FDR reads
			noBlacklist = bedIn.replace('.bed', '_noBlacklist.bed')
			# Intersect input reads with the blacklist region, and return those that do not intersect 
			outfh = open(noBlacklist, 'w')
			proc = subprocess.Popen([program, '-a', bedIn, '-b', blacklistregions, '-v'],stdout=outfh)
			proc.communicate()
			outfh.close()
			blackListed=blackListed+[noBlacklist]
		return (blackListed)
	except:
		print "Problem with blacklist."

def runRepeatMask(mappedReads):
	# Usage: Remove repeat regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools) and blastlist regions removed.
	# Output: Bedfile with repeat regions removed.
	program = 'intersectBed'
	repeatregions = os.getcwd() + '/docs/repeat_masker.bed'
	masked=[]
	try:
		for bedIn in mappedReads:
			# File name for low FDR reads
			noRepeat = bedIn.replace('.bed', '_noRepeat.bed')
			# Intersect input reads with the repeat region, and return those that do not intersect 
			outfh = open(noRepeat, 'w')
			proc = subprocess.Popen([program, '-a', bedIn, '-b', repeatregions, '-v'],stdout=outfh)
			proc.communicate()
			outfh.close()
			masked=masked+[noRepeat]
		return (masked)
	except:
		print "Problem with repeat masking."

def extractExons(bedIn):
	# Usage: Extract all exonic reads from a bedfile
	# Input: .bed file
	# Output: Exonic bedfile
	program = 'intersectBed'
	exonsBed = os.getcwd() + '/docs/allExons.bed'
	try:
		# File name for low FDR reads
		exonicReads = bedIn.replace('.bed', '_exons.bed')
		# Intersect input reads with the blacklist region, and return those that do not intersect 
		outfh = open(exonicReads, 'w')
		proc = subprocess.Popen([program, '-a', bedIn, '-b', exonsBed],stdout=outfh)
		proc.communicate()
		outfh.close()
		return (exonicReads)
	except:
		print "Problem with extraction of exonic reads."

def extractIntrons(bedIn):
	# Usage: Extract all intronic reads from a bedfile
	# Input: .bed file
	# Output: Intronic bedfile and all non-intronic reads
	program = 'intersectBed'
	intronsBed = os.getcwd() + '/docs/allIntrons.bed'
	try:
		# File name for low FDR reads
		intronicReads = bedIn.replace('.bed', '_introns.bed')
		nonIntronicReads = bedIn.replace('.bed', '_Notintrons.bed')
		# Intersect input reads with the intronic bedFile and return those that intersect 
		outfh = open(intronicReads, 'w')
		# proc = subprocess.Popen([program, '-a', bedIn, '-b', intronsBed, '-wa'],stdout=outfh)
		proc = subprocess.Popen([program, '-a', bedIn, '-b', intronsBed,'-u'],stdout=outfh)
		proc.communicate()
		outfh.close()
		# Return all reads that do not intersect
		outfh = open(nonIntronicReads, 'w')
		# proc = subprocess.Popen([program,'-a', bedIn,'-b', intronsBed,'-wa','-v'],stdout=outfh)
		proc = subprocess.Popen([program,'-a', bedIn,'-b', intronsBed,'-v'],stdout=outfh)
		proc.communicate()
		outfh.close()
		# Return file handles 
		return (intronicReads,nonIntronicReads)
	except:
		print "Problem with extraction of intronic reads."

def extractUTRs(bedIn):
	# Usage: Extract all UTR specific reads from the input file.
	# Input: .bed file
	# Output: Mutually exclusive partitions of the input file.
	program = 'intersectBed'
	fivePUTRBed = os.getcwd() + '/docs/5pUTRs_Ensbl_sort_clean_uniq.bed'
	threePUTRBed = os.getcwd() + '/docs/3pUTRs_Ensbl_sort_clean_uniq.bed'
	cdsBed = os.getcwd() + '/docs/Exons_Ensbl_sort_clean_uniq.bed'
	try:
		# Extract 5p reads and NOT 5p reads 
		fivePreads = bedIn.replace('.bed', '_5p.bed')
		notFivePreads = bedIn.replace('.bed', '_NOT5p.bed')
		outfh = open(fivePreads, 'w')
		proc = subprocess.Popen([program, '-a', bedIn, '-b', fivePUTRBed,'-u','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		outfh = open(notFivePreads, 'w')
		proc = subprocess.Popen([program, '-a', bedIn, '-b', fivePUTRBed,'-v','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()

		# Extract 3p UTR reads and NOT 3pUTR
		threePreads = bedIn.replace('.bed', '_3p.bed')
		notThreePreads = bedIn.replace('.bed', '_NOT3p.bed')
		outfh = open(threePreads, 'w')
		proc = subprocess.Popen([program, '-a', notFivePreads, '-b', threePUTRBed,'-u','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		outfh = open(notThreePreads, 'w')
		proc = subprocess.Popen([program, '-a', notFivePreads, '-b', threePUTRBed,'-v','-s'],stdout=outfh)
		proc.communicate()

		# Extract CDS reads and NOT CDS reads 
		CDSreads = bedIn.replace('.bed', '_cds.bed')
		notCDSreads = bedIn.replace('.bed', '_NOTcds.bed')
		outfh = open(CDSreads, 'w')
		proc = subprocess.Popen([program, '-a', notThreePreads, '-b', cdsBed,'-u','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		outfh = open(notCDSreads, 'w')
		proc = subprocess.Popen([program, '-a', notThreePreads, '-b', cdsBed,'-v','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		
		outfh.close()
		return (fivePreads,notFivePreads,CDSreads,notCDSreads,threePreads,notThreePreads)
	except:
		print "Problem with extraction of UTR reads."

def sortCLIPClusters(CLIPPERclusters):
	# Usage: From a bedfile of CLIPper clusters, sort them based upon read count
	# Input: BedFile
	# Output: Sorted bedFile
	try:
		# Name and open the output file
		sortedClusters = CLIPPERclusters.replace('.bed','_sortedClusters.bed')
		# Open file for reading and writing (if file does not exist, it will create one.)
		outfh = open(sortedClusters, 'w+')
		# Process the output file 
		with open(CLIPPERclusters, 'r') as infile:
			for line in infile:	
				chrom,start,end,name,stats,strand = line.strip().split('\t')
				readCountPerCluster = name.strip().split('_')[2]
				outfh.write('\t' .join((chrom,start,end,name,readCountPerCluster,'\n')))
		outfh.close 

		# Extract all clusters	
		outfh = open(sortedClusters)
		cluster_list = [line.strip().split('\t') for line in outfh]
		# Sort using a the read count field
		cluster_list_sorted = sorted(cluster_list, key=lambda line: (int(line[4])), reverse=True)
		outfh.close 

		# Write to the output file
		with open(sortedClusters,'w') as fout:
			for cluster in cluster_list_sorted:
				fout.write('{0}\n'.format('\t'.join(cluster)))
		fout.close 
		
 		return sortedClusters

	except:
		print "Error sorting CLIPper clusters."

def shuffleBedFile(inBed):
	# Usage: Shuffle a bed file (usually used a background file for HOMER analysis) 
	# Input: Bedfile
	# Output: SHuffled bedFile
	program = 'shuffleBed'
	referenceFile = os.getcwd()+'/docs/hg19_transcriptome_collapse_exon.bed'
	genomeFile = os.getcwd()+'/docs/human.hg19.genome'
	try:
		shuffledBed = inBed.replace('.bed','_shuffled.bed')
		outfh = open(shuffledBed, 'w')
		proc = subprocess.Popen([program,'-i',inBed,'-incl',referenceFile,'-g',genomeFile],stdout=outfh)
		proc.communicate()
		return shuffledBed
	except:
		print "Problem generating shuffled bedfile."

def makeBedForHOMER(inBed):
	# Usage: This modified the bedfile for processing homer by making the first field a concatenation of chr_start_end
	# Input: Clean bedfile with first 5 field properly assigned
	# Output: Modified bed file with first field chr_start_end, and name excluded 
	try:
		# Make sure bedfile only has 5 fields
		bedForHOMER=inBed.replace('.bed','_forHOMER.bed')	
		# Open new file
		f = open(bedForHOMER, 'w')
		with open(inBed, 'r') as infile:
			for line in infile:	
				elementList = line.strip().split('\t')
				# Re-write the bed file with chr replaced
				f.write('\t' .join((elementList[0]+'_'+elementList[1]+'_'+elementList[2],elementList[0],elementList[1],elementList[2],elementList[5],'\n')))
		f.close()
		return bedForHOMER
	except:
		print "Error making bed file for HOMER."

def runHOMER(inBed,outDirName):
	# Usage: Run the HOMER motif finding algorithm 
	# Input: Bedfile properly modified for HOMER
	# Output: A directory containing the HOMER output files
	program='findMotifsGenome.pl'
	program2='annotatePeaks.pl'
	homerReferenceFile = os.getcwd()+'/docs/hg19_transcriptome_collapse_exon.bed'
	try:
		# Convert the input bedFile into HOMER compatible format
		inBedForHOMER=makeBedForHOMER(inBed)
		# Get the path of the input file 
		path,filename=os.path.split(inBedForHOMER)
		outDir=path+'/'+outDirName
		# Call HOMER, which will generate a directory of files
		proc = subprocess.Popen([program,inBedForHOMER,'hg19',outDir,'-rna','-bg',homerReferenceFile])
		proc.communicate()
	except:
		print "Error running HOMER on: %s" %outDirName

def getReadDist(geneList):
	# Usage: Extract raw reads and associated gene start/stop associated with each Ensembl gene ID
	# Input: Gene list
	# Output: Bed file with read IDs and assocaited start/stop for the gene
	geneStartStopRepo = os.getcwd() + '/docs/all_genes.txt'
	# try:
	# Reads associated with each gene type
	geneReadList=geneList+'_LowFDRreads.bed'
	# Ensure that we read genes from the filtered lncRNAs
	if 'lincRNA' in geneList:
		geneReadList=outfilepath+'/clipGenes_lincRNA_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed'
	readsByGenePosition=geneList+'_genePosition'
	# Open the output file
	outfh = open(readsByGenePosition, 'w')
	with open(geneList, 'r') as infile:
		for ensemblID in infile:
			# Gene name
			ensemblID=ensemblID.strip().split('\t')[0]
			# Extract basic gene data (from the first of the list, if there are repeats)
			storePosition=grep(ensemblID,geneStartStopRepo)
			if storePosition:
				# Parse information about the gene
				ID,chrom,geneStart,geneEnd=storePosition[0].strip().split('\t')
				# Then extract all raw reads with the gene ID and output data for each read
				results=grep(ensemblID,geneReadList)
				# Check to make sure reads were recovered
				if results:
					for read in results:
						readData=read.strip().split('\t')
						outfh.write('\t'.join((ID,geneStart,geneEnd,readData[1],readData[2],'\n')))
	return readsByGenePosition
	# except:
	# 	print "Error extracting read position per gene."

def initializeParameters():
	# < THIS WILL BE USED IN LATER VERSION OF THE SCRIPT > # 
	parser = OptionParser()
	parser.add_option("--proj", dest="project",help="Name of input project.")
	parser.add_option("--genome", dest="genome",help="Bowtie genome index.")	
	parser.add_option("--seed", dest="seedLine",help="Line in seedfile to read.")
	# Initialize option parser
	(options,args) =  parser.parse_args()
	# Error message
	error = ""
	# Handle the exceptions
	if not options.project:
		error += "\nError: Please specify a CLIP project."
	if not options.seedLine:
		error += "\nError: Please specify a line in seedfile to read."
	if options.genome:
		genomeIndex=str(options.genome)
	print error
	# Return the name of the project
	return (str(options.project),int(options.seedLine))

def filterSnoRNAs(proteinCodingReads):
	# Usage: Filter snoRNA and miRNAs from protein coding reads.
	# Input: .bed file with protein coding reads.
	# Output: snoRNA and miR filtered .bed file.
	program = 'intersectBed'
	snoRNAindex = os.getcwd() + '/docs/snoRNA_reference/snoRNAmasker_formatted.bed'
	miRNAindex = os.getcwd() + '/docs/miR_sort_clean.bed'

	try:
		# File name for low FDR reads
		proteinWithoutsnoRNAs = proteinCodingReads.replace('.bed', '_snoRNAremoved.bed')
		proteinWithoutmiRNAs = proteinWithoutsnoRNAs.replace('.bed', '_miRNAremoved.bed')

		# Intersect protein coding reads with snoRNA masker
		outfh = open(proteinWithoutsnoRNAs, 'w')
		proc = subprocess.Popen([program,'-a', proteinCodingReads,'-b', snoRNAindex,'-v','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		# Intersect protein coding reads with miR masker
		outfh = open(proteinWithoutmiRNAs, 'w')
		proc = subprocess.Popen([program,'-a', proteinWithoutsnoRNAs,'-b', miRNAindex,'-v','-s'],stdout=outfh)
		proc.communicate()
		outfh.close()
		return (proteinWithoutmiRNAs)
	except:
		logOpen.write("Problem obtaining lowFDR reads.\n")
		print "Problem obtaining lowFDR reads."

def convertENBLids(inNames):
	# Usage: Converl ENST to ENSG (unique ID) using ENSEMBL annotation file
	# Input: List of ENST IDs
	# Output: List of ENSG IDs

	genesFile = os.getcwd() + '/docs/hg19_ensembl_genes.txt'
	ensemblIDfile=np.genfromtxt(genesFile,usecols=(1,12,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
	temp=[]
	for name in inNames:
		outName=ensemblIDfile[ensemblIDfile[:,0]==name,1]
		temp=temp+[outName]
	temp=np.array(temp)
	return temp

def parseRepeatMapped(RTmerged):
	# Usage: Parse repeat bed file for each 
	# In: bedFile with merged RT stops
	# Out: List of lists (each has RT stops assocaited with a single gene)
	
	# Load start coordinates
	start=np.loadtxt(RTmerged,usecols=(1,),delimiter='\t',dtype=int)
	start[start < 0]=0 # Ensure no negative start coordinates and normalize 
	rtStops=start+15 # Regenerate the RT stop position: RTmerged is +/- 15 bases from the RT stop
	rnu_u1=rtStops[rtStops < 165]
	rnu_u2=rtStops[(rtStops > 164) & (rtStops < 353)]
	rnu_u4_1=rtStops[(rtStops > 352) & (rtStops < 497)]
	rnu_u4_atac=rtStops[(rtStops > 496) & (rtStops < 627)]
	rnu_u5=rtStops[(rtStops > 626) & (rtStops < 744)]
	rnu_u6=rtStops[(rtStops > 743) & (rtStops < 850)]
	rnu_u6_atac=rtStops[(rtStops > 849) & (rtStops < 975)]
	rnu_u7=rtStops[(rtStops > 974) & (rtStops < 1038)]
	rnu_u11=rtStops[(rtStops > 1037) & (rtStops < 1173)]
	rnu_u12=rtStops[(rtStops > 1172) & (rtStops < 1322)]
	rn7_sl2=rtStops[(rtStops > 1321) & (rtStops < 1621)]
	rn7_sk=rtStops[(rtStops > 1620) & (rtStops < 1953)]
	rny_1=rtStops[(rtStops > 1952) & (rtStops < 2066)]
	rny_3=rtStops[(rtStops > 2065) & (rtStops < 2168)]
	rny_4=rtStops[(rtStops > 2167) & (rtStops < 2264)]
	rny_5=rtStops[(rtStops > 2263) & (rtStops < 2348)]
	u3a=rtStops[(rtStops > 2347) & (rtStops < 2565)]
	rna_5s=rtStops[(rtStops > 2564) & (rtStops < 2686)]
	
	# Define the end of the rRNA regions
	global rRNAstart
	rRNAstart=2686
	rDNA=rtStops[(rtStops > rRNAstart)]
	readLists=[rnu_u1,rnu_u2,rnu_u4_1,rnu_u4_atac,rnu_u5,rnu_u6,rnu_u6_atac,rnu_u7,rnu_u11,rnu_u12,rn7_sl2,rn7_sk,rny_1,rny_3,rny_4,rny_5,u3a,rna_5s,rDNA]
	
	# Save the files 
	labels=['U1','U2','U4_1','U4_ATAC','U5','U6','U6_ATAC','U7','U11','U12','7SL-SRP','7SK','Y_1','Y_3','Y_4','Y_5','U3A','5s','rDNA']
	for i in range(len(readLists)):
		repeatData=outfilepath+'SourceData_RepeatRNA_%s'%labels[i]
		np.savetxt(repeatData,readLists[i],delimiter='\t',fmt="%s")

	return readLists

def makeGeneLists(outfilepath,CLIPPERlowFDR):
	# Usage: Make sorted gene lists for each RNA type, using specific reference file for snoRNAs and filtering miRNA, snoRNA from protein coding and lincRNAs.
	# Input: Outfile path and low FDR reads.
	# Output: Path to each file, number of reads assigned to each gene type, list of labels for the order gene types analyzed.

	filesToSave=[]
	readsPerGeneType=[]
	geneTypeLabels=[]
	# Get all lowFDR files by each gene type
	for bedFile in glob.glob(outfilepath+"*_LowFDRreads.bed"):
		head, tail = os.path.split(bedFile)
		geneType=tail.strip().split('_')[1]
		print geneType
		geneTypeLabels=geneTypeLabels+[geneType]
		# For protein coding genes ...
		if geneType == "proteinCoding":
			# ... we will use the center coordinate file with snoRNAs filtered
			bedFile=outfilepath+"clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed"
			geneNames=np.genfromtxt(bedFile,usecols=(3,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
			nameList=np.array([i.strip().split('.')[0] for i in geneNames],dtype='string')
		# For lincRNAs genes (in which snoRNAs are contained) ...
		elif geneType == "lincRNA":
			# ... we will use the center coordinate file with snoRNAs filtered
			bedFile=outfilepath+"clipGenes_lincRNA_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed"
			geneNames=np.genfromtxt(bedFile,usecols=(3,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
			nameList=np.array([i.strip().split('.')[0] for i in geneNames],dtype='string')
		# For snoRNAs ...
		elif geneType == "snoRNA":
			# ... we will intersect centerpoints of all lowFDR reads with a reference file of well-annotated (accepted) snoRNAs
			program = 'intersectBed'
			snoRNAindex =  os.getcwd() + '/docs/snoRNA_reference/sno_coordinates_hg19_formatted.bed'
			# This has the correct 6 column bedFile format			
			bedFile = outfilepath+'clipGenes_snoRNA_LowFDRreads.bed'
			outfh = open(bedFile, 'w')
			# This produces the primary bedGraph / bigWig for visualization
			CLIPPERlowFDRcenters=getBedCenterPoints(CLIPPERlowFDR)
			allLowFDRCentersBedGraph=makeBedGraph(CLIPPERlowFDRcenters)
			# -u reports only one overlap for each read in -a, and -s reports strand
			proc = subprocess.Popen([program,'-a', CLIPPERlowFDRcenters, '-b', snoRNAindex,'-wb'],stdout=outfh)
			proc.communicate()
			outfh.close()	
			snoRNAbedGraph=makeBedGraph(bedFile)
			# Unique snoRNA ID is in column of the resulting file 
			nameList=np.genfromtxt(bedFile,usecols=(9,),delimiter='\t',dtype='string')
			# nameList=np.array([i.strip().split('.')[0] for i in geneNames],dtype='string')
		else:
			# All other lowFDR read files (extracted using Ensembl ID) have gene name in the ninth column 
			geneNames=np.genfromtxt(bedFile,usecols=(9,),delimiter='\t',dtype='string')
			nameList=np.array([i.strip().split('.')[0] for i in geneNames],dtype='string')
		
		# Count number of reads per gene name 
		try:	
			# Count the number of occourances of each gene name
			count={}
			for name in nameList:
				if count.has_key(name):
					count[name] += 1
				else:
					count[name] = 1
			# Sort the dictionary by value
			countSort=sorted(count.items(),key=lambda item: item[1])
		except:
			print "No genes of this type."
		
		# Save the sorted gene names
		tosave=outfilepath+'clipGenes_%s'%geneType
		outfh = open(tosave, 'w')
	 	for pair in countSort:
	 		string=pair[0]+'\t'+str(pair[1])+'\n'
	 		outfh.write(string)
	 	outfh.close()

	 	# Store files
	 	filesToSave=filesToSave+[tosave]

	 	# Save list using "SourceData" handle for easy extraction 
	 	tosave=outfilepath+'SourceData_geneList_%s'%geneType
		outfh = open(tosave, 'w')
	 	for pair in countSort:
	 		string=pair[0]+'\t'+str(pair[1])+'\n'
	 		outfh.write(string)
	 	outfh.close()

	return(filesToSave)

def getBedCenterPoints(inBed):
	# Usage: Obtain ceter coordiantes of bedFile
	# Input: BedFile.
	# Output: Center coodinates returned.
	
	try:
		# Make sure bedfile only has 5 fields
		outBed=inBed.replace('.bed','_centerCoord.bed')	
		# Open new files
		f = open(outBed, 'w')
		with open(inBed, 'r') as infile:
			for line in infile:	
				elementList = line.strip().split('\t')
				# Re-write the CLIPper windows file
				f.write('\t' .join((elementList[0],str(int(elementList[1])+15),str(int(elementList[1])+16),elementList[9],elementList[4],elementList[5],'\n')))
		f.close()
		return outBed
	except:
		print "Problem with center coordiante extacrtion."

def partitionSnoRNAs(outfilepath):
	# Usage: Plot abundance of reads that map to each snoRNA type and coverage histograms for each.
	# Input: Output file path.
	# Output: Figure.

	program='intersectBed'
	snoRNAindex = os.getcwd() + '/docs/snoRNA_reference/sno_coordinates_hg19_formatted.bed'
	# File name for low FDR snoRNA reads
	snoRNAreads = outfilepath+'clipGenes_snoRNA_LowFDRreads.bed'
	# snoRNA read types 
	hcBOX_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_HCAbox.bed'
	SCA_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_SCARNA.bed'
	CDBOX_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_CDBOX.bed'
	outfh1 = open(hcBOX_READS, 'w')
	outfh2 = open(SCA_READS, 'w')
	outfh3 = open(CDBOX_READS,'w')

	# Read file into independent files for each snoRNA type
	scaRNAs=0
	hcaBox=0
	cdBox=0
	with open(snoRNAreads, 'r') as infile:
		for line in infile:	
			elementList = line.strip().split('\t')
			snoType = elementList[10]
			elementList.append('\n')
			if snoType == 'C':
				outfh3.write('\t'.join(elementList))
				cdBox += 1
			elif snoType == 'H':
				outfh1.write('\t'.join(elementList))
				hcaBox += 1
			elif snoType == 's':
				outfh2.write('\t'.join(elementList))
				scaRNAs += 1
	# Close the input files before moving forward and counting reads
	outfh1.close()
	outfh2.close()
	outfh3.close()

	# Plot pie chart 
	fig6=plt.figure(6)
	ax=plt.subplot(2,2,1)
	totalCount=float(scaRNAs+hcaBox+cdBox)
	
	# Check to make sure snoRNA reads recovered
	if totalCount > 0:

		# Make a dictionary of gene start,stop for each unique snoRNA ID
		d={}
		snoStartStopRepo =  os.getcwd() + '/docs/snoRNA_reference/sno_coordinates_hg19_formatted.bed'
		with open(snoStartStopRepo, 'r') as infile:
			for gene in infile:	
				chrom,geneStart,geneEnd,name,snoType,strand=gene.strip().split('\t')
				d[name]=[geneStart,geneEnd]
			infile.close()
		
		# Read each snoRNA file
		snoRNAfiles=[CDBOX_READS,hcBOX_READS,SCA_READS]
		# Store for position files
		positionFiles=[]
		for snoFile in snoRNAfiles:
			# Create a list of the output files
			readsByGenePosition=snoFile.replace('.bed', '_genePosition_snoRNA')
			positionFiles=positionFiles+[readsByGenePosition]
			# Open the output file
			outfh = open(readsByGenePosition, 'w')
			# Read through the snoRNA file
			count=0
			with open(snoFile, 'r') as infile:
				for read in infile:
					# Extract the elements of each read
					elementList=read.strip().split('\t')
					# Extract the snoRNAID
					snoID=elementList[9].strip()
					# Extract information for each gene
					geneStartStop=d[snoID]
					# Write to the output file: start and end of snoRNA based upon reference file, and start and end of snoRNA read
					outfh.write('\t'.join((snoID,geneStartStop[0],geneStartStop[1],elementList[1],elementList[2],'\n')))
					count += 1
			infile.close()
			outfh.close()

def saveSourceData(UTRgenelists):
	# Usage: Save source data for UTR analysis 
	# In: Path to lists of reads partitioned by type
	# Out: None

	UTRlabelsTosave=['5p','cds','3p']
	j=0
	for UTRgenelist in UTRgenelists:
		# Read all names
		geneNames=np.genfromtxt(UTRgenelist,usecols=(3,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
		nameList=np.array([i.strip().split('.')[0] for i in geneNames],dtype='string')

		# Count reads per name
		count={}
		for name in nameList:
			if count.has_key(name):
				count[name] += 1
			else:
				count[name] = 1
		# Sort the dictionary by value
		countSort=sorted(count.items(),key=lambda item: item[1])

		# Save the sorted gene names
		tosave=outfilepath+'SourceData_ReadsByUTR_%s' %UTRlabelsTosave[j]
		outfh = open(tosave, 'w')
		for pair in countSort:

		 	string=pair[0]+'\t'+str(pair[1])+'\n'
		 	outfh.write(string)
		outfh.close()
		j += 1

def plot_ReadAccounting(outfilepath,sampleName):
	# Usage: Make hbar chart of read count at specific steps in pipeline
	# Input: File path and sample name
	# Output: None

	# Paths to files for which we want to count reads
	rawRead1=outfilepath+sampleName+'_R1.fastq'
	rawRead2=outfilepath+sampleName+'_R2.fastq'
	reads3pTrim=[outfilepath+sampleName+'_R1_3ptrimmed.fastq',outfilepath+sampleName+'_R2_3ptrimmed.fastq']
	readsFilter=[outfilepath+sampleName+'_R1_3ptrimmed_filter.fastq',outfilepath+sampleName+'_R2_3ptrimmed_filter.fastq']
	readsNoDupes=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe.fastq',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe.fastq']
	readsMappedReapeat=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_mappedToRepeatRNA.bed',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_mappedToRepeatRNA.bed']
	readsMappedHg19=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes.bed',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes.bed']
	readsMappedBlacklist=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes.bed',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes.bed']
	readsMappedRepeatMask=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes_noBlacklist_noRepeat.bed',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedToRepeat_mapped_nodupes_noBlacklist_noRepeat.bed']
	clipperIN=outfilepath+sampleName+'_allreads.mergedRT_CLIPPERin.bed'
	clipperOUT=outfilepath+sampleName+'_allreads.mergedRT_CLIP_clusters_lowFDRreads.bed'
	# File name labels for plotting
	fileNames=['Raw (R1)','Raw (R2)','3p Trim (R1)','3p Trim (R2)','Filter (R1)','Filter (R2)','No dupes (R1)','No dupes (R2)','RepeatMapped(R1)','RepeatMaped(R2)','Hg19Mapped (R1)','Hg19Mapped(R2)','Blacklist (R1)','Blacklist (R2)','RepeatMask(R1)','RepeatMask(R2)','ClipperIn','ClipperOut']
	# Determine read count for specified files
	filesToCount=[rawRead1,rawRead2,reads3pTrim[0],reads3pTrim[1],readsFilter[0],readsFilter[1],readsNoDupes[0],readsNoDupes[1],readsMappedReapeat[0],readsMappedReapeat[1],readsMappedHg19[0],readsMappedHg19[1],readsMappedBlacklist[0],readsMappedBlacklist[1],readsMappedRepeatMask[0],readsMappedRepeatMask[1],clipperIN,clipperOUT]
	counts=[]
	counter=0
	for fileString in filesToCount:
		temp=lineCount(fileString)
		if counter < 8:
			temp=temp/4  # Ensure that the read count for all .fastq files is corrected (1 read = 4 lines)
		counts=counts+[temp]
		counter += 1
	# Plot
	#fig1=plt.figure(1)
	plt.subplot(2,3,1)
	ind = np.arange(len(counts)) + 0.5
	plt.barh(ind,list(reversed(counts)),align='center',color='blue')
	plt.xlabel('Read count per file',fontsize=5)
	plt.tick_params(axis='xticks',labelsize=5) 
	plt.yticks(ind,list(reversed(fileNames)),fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5) 
	# Remove tick lines 
	ax=plt.gca()
	for line in ax.get_yticklines():
	    line.set_markersize(0)
	plt.title('Read counts per step',fontsize=10)
	# Save plot data 
	readDF=pd.DataFrame()
	readDF['File_name']=fileNames
	readDF['Reads_per_file']=counts
	outfilepath=outfilepath + 'SourceData_plot_ReadAccounting'
	readDF.to_csv(outfilepath)

def plot_BoundGeneTypes(outfilepath,sampleName):
	# Usage: Quantify bound genes of select type
	# Input: Gene lists to plot, file path
	# Output: None

	# Gene lists
	geneListToPlot=[outfilepath+'clipGenes_proteinCoding',outfilepath+'clipGenes_pseudogenes',outfilepath+'clipGenes_snoRNA',outfilepath+'clipGenes_lincRNA',outfilepath+'clipGenes_miR']
	# Count the number of genes in each list
	counts=[]
	labels=[]
	for fileString in geneListToPlot:
		head, tail = os.path.split(fileString)
		geneTypeName = tail.split("_")[1]
		labels=labels+[geneTypeName]
		temp=lineCount(fileString)
		counts=counts+[temp]
	# Sort the counts 
	labels=np.array(labels)
	counts=np.array(counts)
	sortedIndex=[i[0] for i in sorted(enumerate(counts),key=lambda x:x[1])]
	# Plot
	#fig1=plt.figure(1)
	plt.subplot(2,3,6)
	ind = np.arange(len(counts)) + 0.5
	plt.bar(ind,counts[sortedIndex],align='center',color='blue')
	locs,pltlabels = plt.xticks(ind,labels[sortedIndex],fontsize=5)
	plt.setp(pltlabels, rotation=90, fontsize=5)
	plt.tick_params(axis='xticks',labelsize=5) 
	ax=plt.gca()
	for line in ax.get_xticklines():
		line.set_markersize(0)
	plt.ylabel('Gene count',fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5)
	plt.title('Read counts per step',fontsize=10)
	plt.title('Bound gene type and count',fontsize=10)
	# Save plot data 
	readDF=pd.DataFrame()
	readDF['Gene_type']=labels[sortedIndex]
	readDF['Number_bound_genes']=counts[sortedIndex]
	outfilepath=outfilepath + 'SourceData_plot_BoundGeneTypes'
	readDF.to_csv(outfilepath)

def plot_ReadsPerCluster(outfilepath,sampleName):
	# Usage: Make histogram of RT stops per cluster
	# Input: File path and sample name
	# Output: None

	# Extract file containing the number of reads per cluster
	CLIPpeReadsPerCluster=outfilepath+sampleName+'_allreads.mergedRT_CLIP_clusters.readsPerCluster'
	CLIPpeReadsPerClusterList=[]
	with open(CLIPpeReadsPerCluster) as fin:
		CLIPpeReadsPerClusterList = [int(line.strip()) for line in fin]
	# Make hisogram, with each bin (except last) including lower value, excluding upper
	bins=range(min(CLIPpeReadsPerClusterList)-10,max(CLIPpeReadsPerClusterList)+10,10)
	hist,bins=np.histogram(CLIPpeReadsPerClusterList,bins=bins)
	# Plot
	#fig1=plt.figure(1)
	plt.subplot(2,3,2) 
	xlimit=100 # The maximum number of reads per cluster visualized 
	width=0.7*(bins[1]-bins[0]) # Bar width is 70% of bin size
	center=(bins[:-1] + bins[1:])/2 # The center position of each bin
	plt.bar(center, hist,align='center',width=width)
	plt.xlabel('Reads per cluster (bin=10)',fontsize=5)
	plt.ylabel('Frequency (RT stop count)',fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5) 
	plt.xlim(0,xlimit)
	plt.title('Reads per cluster',fontsize=10)

def plot_ClusterSizes(outfilepath,sampleName):
	# Usage: Make boxplot of cluster lengths
	# Input: File path
	# Output: None

	# Cluster file returned from clipper
	clipClusters=outfilepath+sampleName+"_allreads.mergedRT_CLIP_clusters"
	# Obtain lengths of each cluster 
	clusterLengths=[]
	with open(clipClusters) as fin:
		for line in fin:
			elements=line.strip().split('\t')
			start=int(elements[1])
			end=int(elements[2])
			lengthOfCluster=math.fabs(start-end)
			clusterLengths=clusterLengths+[lengthOfCluster]
	# Plot
	#fig1=plt.figure(1)
	plt.subplot(2,3,3)
	ylimit=500 # Maximum cluster length to show 
	plt.boxplot(clusterLengths)
	plt.tick_params(axis='x',labelbottom='off') 
	ax=plt.gca()
	for line in ax.get_xticklines():
		line.set_markersize(0)
	plt.ylabel('Cluster length (bases)',fontsize=5)
	plt.tick_params(axis='yticks',labelsize=5)
	plt.ylim(0,ylimit)
	plt.title('Cluster size',fontsize=10)

def plot_clusterBindingIntensity(outfilepath,sampleName):
	# Usage: Read distribution around cluster center point
	# Input: File path
	# Output: None 

	# Get heatmap of read intensity around cluster centers	
	clusterCenterHeatmap=outfilepath+sampleName+'_allreads.mergedRT_CLIP_clusters_cleaned_sorted.clusterCenter_heatmap.txt'
	heatmapCols=40
	# Skip the header row (with column numbers) and the first column (cluster coordinates)
	heatmapData=np.loadtxt(clusterCenterHeatmap,skiprows=1,dtype='float',usecols=range(1,heatmapCols+1))
	# Compute the sum of each row, where each row is centered on the coordiantes of each cluster center
	clusterSums=heatmapData.sum(axis=1)
	# Sort by sum of read abundance per cluster
	sortedIndex=[i[0] for i in sorted(enumerate(clusterSums),key=lambda x:x[1])]
	sortedHeatMap=heatmapData[sortedIndex,:]
	# Plot
	#fig1=plt.figure(1)
	plt.subplot(2,3,4)
	ylimit=sortedHeatMap.shape[0]
	# The final value in the matrix (largest, in this case) is plotted at the top of the heatmap
	p=plt.pcolormesh(sortedHeatMap,cmap='Blues')
	plt.colorbar(p)
	plt.xticks(range(sortedHeatMap.shape[1]))
	plt.tick_params(axis='x',labelbottom='off') 
	plt.xlabel('Cluster position',fontsize=5)
	plt.ylim(0,ylimit)
	plt.ylabel('Cluster number',fontsize=5)
	plt.title('Read distribution',fontsize=10)

def plot_readsBymRNAregion(outfilepath,sampleName):
	utrData=['5p','cds','3p'] 
	readCounts=[]
	for utr in utrData:
		# Get read data from source files (sorted gene list) and check against raw .bed file
		source=outfilepath + 'SourceData_ReadsByUTR_%s'%utr
		geneCount=np.genfromtxt(source,usecols=(1,),delimiter='\t',dtype='int')
		bedSource=outfilepath + 'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved_%s.bed'%utr
		bedFile=np.genfromtxt(bedSource,usecols=(3,),delimiter='\t',dtype='string')
		readCounts=readCounts+[bedFile.size]
	# Plot
	# Pie chart based upon read based accounting
	allProteinCoding=outfilepath + 'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed'
	allProteinCoding=np.genfromtxt(allProteinCoding,usecols=(3,),delimiter='\t',dtype='string')
	totalReads=allProteinCoding.size
	nonIntronic=sum(readCounts) # Reads that fall into 5p UTR, CDS, 3p UTR
	vector=[float(totalReads-nonIntronic)/totalReads,float(readCounts[0])/totalReads,float(readCounts[1])/totalReads,float(readCounts[2])/totalReads]
	ax=plt.subplot(2,3,5)
	pie_wedges=ax.pie(vector,labels=["Intronic","5p UTR","CDS","3pUTR"],labeldistance=1.1,autopct='%1.1f%%')
	plt.rcParams['font.size']=10
	for wedge in pie_wedges[0]:
		wedge.set_edgecolor('black')
		wedge.set_lw(1)

def plot_mRNAgeneBodyDist(outfilepath,sampleName):
	# Bed file with protein coding reads
	filteredProteinCoding = outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed'
	# Generate gene-by-gene read coverage histogram 
	averageGraph=outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved_cleaned_sorted_UTRs_scaled_cds200_abt0_averageGraph.txt'
	# Number of columns 
	avgGraphCols=600
	avgGraphData=np.loadtxt(averageGraph,skiprows=1,dtype='float',usecols=range(1,avgGraphCols+1))
	fig2=plt.figure(2)
	plt.subplot2grid((2,3),(0,0),colspan=3)
	ylimit=max(avgGraphData[1,:])*1.1
	plt.plot(avgGraphData[1,:],color='blue',linewidth='2')
	plt.ylim(0,ylimit)
	plt.vlines(200,0,ylimit,linestyles='dashed')
	plt.vlines(400,0,ylimit,linestyles='dashed')
	plt.tick_params(axis='x',labelbottom='off') 
	plt.xlabel('mRNA gene body (5pUTR, CDS, 3pUTR)')
	plt.ylabel('Abundance')
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('CLIP signal across average mRNA transcript.',fontsize=10)

def convertENBLids(inNames):
	# Usage: Converl ENST to ENSG (unique ID) using ENSEMBL annotation file
	# Input: List of ENST IDs
	# Output: List of ENSG IDs
	genesFile = os.getcwd() + '/docs/hg19_ensembl_genes.txt'
	ensemblIDfile=np.genfromtxt(genesFile,usecols=(1,12,),delimiter='\t',dtype='string') # Note that column lookup is zero indexed
	temp=[]
	for name in inNames:
		outName=ensemblIDfile[ensemblIDfile[:,0]==name,1]
		temp=temp+[outName]
	temp=np.array(temp)
	return temp

def extractGenesByUTR(utrType,treatMatrixData,geneNames):
	# Usage: Extract genes with exclusive binding 
	# Input: UTR type, data, names
	# Output: Extracted and sorted data, names
	# Extract specific regions
	if utrType=='5p':
		indexer=(treatMatrixData[:,range(200,600)].sum(axis=1) == 0) & (treatMatrixData[:,range(0,200)].sum(axis=1) > 0)
	elif utrType=='3p':
		indexer=(treatMatrixData[:,range(0,400)].sum(axis=1) == 0) & (treatMatrixData[:,range(400,600)].sum(axis=1) > 0)
	else:
		indexer=(treatMatrixData[:,range(0,200)].sum(axis=1) == 0) & (treatMatrixData[:,range(400,600)].sum(axis=1) == 0) & (treatMatrixData[:,range(200,400)].sum(axis=1) > 0)
	# Extract data    
	extractedData=treatMatrixData[indexer,:]
	extractedNames=geneNames[indexer]
	# Get total abundance / gene
	geneSums=extractedData.sum(axis=1)  
	extractedDataSort=extractedData[np.argsort(geneSums),:]
	extractedNamesSort=extractedNames[np.argsort(geneSums)]
	return (extractedDataSort,extractedNamesSort)

def plot_geneBodyPartition(outfilepath,sampleName):
	# Usage: 
	# In: Outfilepath
	# Out: Source data, heatmaps

	# Binding across gene body (5p UTR, CDS, 3p UTR) of protein coding genes
	treatMatrixCols=600
	treatMatrix=outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved_cleaned_sorted_UTRs_scaled_cds200_abt0_treatmatrix.txt'
	treatMatrixData=np.genfromtxt(treatMatrix,skip_header=1,usecols=range(1,treatMatrixCols+1),delimiter='\t',dtype='float')
	geneNames=np.loadtxt(treatMatrix,dtype='string',skiprows=1,usecols=(0,),delimiter='\t')
	# Convert to ENSG IDs
	geneNames=convertENBLids(geneNames)
	# Ensure that genes are in the CLIPper protein coding list (treat analysis can return some spurious genes)
	masterList = outfilepath+'clipGenes_proteinCoding'
	masterNames = np.genfromtxt(masterList,usecols=(0,),delimiter='\t',dtype='string')
	indexer=[]
	for geneName in geneNames:
		if geneName in masterNames:
			indexer=indexer+[1]
		else:
			indexer=indexer+[0]
	indexer=np.array(indexer,dtype=bool)
	geneNames=geneNames[indexer]
	treatMatrixData=treatMatrixData[indexer,:]
	# Isolate genes and treatmatrix data by region
	try:
		data_5p,names_5p=extractGenesByUTR('5p',treatMatrixData,geneNames)
	except:
		print "No 5' exclusive reads."
		data_5p=[]
		names_5p=[]
	try:
		data_cds,names_cds=extractGenesByUTR('cds',treatMatrixData,geneNames)
	except:
		print "No cds exclusive reads."
		data_cds=[]
		names_cds=[]
	try:   
		data_3p,names_3p=extractGenesByUTR('3p',treatMatrixData,geneNames)
	except:
		print "No CDS' exclusive reads."
		data_3p=[]
		names_3p=[]

	# All protein coding genes
	masterList = outfilepath+'clipGenes_proteinCoding'
	masterNames = np.genfromtxt(masterList,usecols=(0,),delimiter='\t',dtype='string')
	# All 5p UTR, CDS, 3p UTR genes
	allTreatMatrixGenes=list(geneNames[:,0])
	# Intronic gene exctraction
	tosave=outfilepath + 'SourceData_geneBodyPartition_ExclusiveIntronic' 
	intronicBoundGenes=list(set(masterNames)-set(allTreatMatrixGenes))
	np.savetxt(tosave,np.array(intronicBoundGenes),fmt="%s")

	# Check relative to previously identified genes
	utrBinding_Names_Fig1=['5p','cds','3p'] 
	utrBinding_Names_Fig2=[names_5p[:,0],names_cds[:,0],names_3p[:,0]]
	utrBinding_Data_Fig2=[data_5p,data_cds,data_3p]
	diffs=[]

	for i in range(0,3):
		# Identified genes with 5p UTR, CDS, 3p UTR binding from read-based paritioning
		utr=utrBinding_Names_Fig1[i]
		source=outfilepath + 'SourceData_ReadsByUTR_%s'%utr
		temp1e=np.genfromtxt(source,usecols=(0,),delimiter='\t',dtype=str)
		# Identified genes with 5p UTR, CDS, 3p UTR from gene body analysis
		temp2b=np.unique(utrBinding_Names_Fig2[i]) 
		# Non-overlaps are due to reads that are ambiguously mapped into different 5p UTR, CDS, 3p UTR region on isoforms
		diff=list(set(temp2b)-set(temp1e))
		diffs.append(diff)

	namesFilter=[]
	dataFilter=[]
	geneTypes=['5p','cds','3p'] 
	titles=['5pUTR','cds','3pUTR'] 
	depth=50
	for i in range(0,3):
		# Names of bound genes
		names=utrBinding_Names_Fig2[i]
		data=utrBinding_Data_Fig2[i]
		# Indexer
		indexer=np.array([geneName not in diffs[i] for geneName in names])
		tempNames=names[indexer]
		namesFilter.append(tempNames)
		tempData=data[indexer,:]
		dataFilter.append(tempData)
		# Save data
		tosave=outfilepath + 'SourceData_geneBodyPartition_Exclusive%s'%titles[i] 
		np.savetxt(tosave,np.array(tempData),fmt="%s")
		# Plot
		plt.subplot2grid((2,3),(1,i),colspan=1)
		p=plt.pcolormesh(dataFilter[i][-depth:-1,:],cmap='Blues')
		plt.title(titles[i],fontsize=10)
		plt.vlines(200,0,depth,linestyles='dashed')
		plt.vlines(400,0,depth,linestyles='dashed')
		plt.tick_params(axis='x',labelbottom='off') 
		plt.tick_params(axis='y',labelleft='off') 
		plt.ylim(0,depth)
		plt.ylabel('Ranked genes (highest on bottom)',fontsize=5)
		plt.xticks(visible=False)
		plt.yticks(visible=False)
		plt.title('%s specific genes: %s'%(titles[i],np.unique(namesFilter[i]).size),fontsize=7.5)

def plot_repeatRNA(outfilepath,sampleName):
	repeatGenomeBuild=os.getcwd() + '/docs/repeat/repeatRNA.fa'
	repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
	repeat_genome_bases=repeat_genome[1]
	repeat_genome_size=len(repeat_genome[1])
	# Repeat index positions
	repeatAnnotation=os.getcwd() + '/docs/repeat/Hs_repeatIndex_positions.txt'
	repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
	repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
	# Python list extraction is not end index inclusive; to extract sequence, use end + 1.
	repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 

	# Get all repeat RNA files
	repeatRNA_RTstops=glob.glob(outfilepath+'SourceData_RepeatRNA_*')
	plotDim=math.ceil(math.sqrt(len(repeatRNA_RTstops)))

	for i in range(len(repeatRNA_RTstops)):
		# RT stop file
		filePath=repeatRNA_RTstops[i]
		# Extract RNA ID
		RNAid=filePath.split('RepeatRNA_')[1]
		# Get start,end for repeat RNA and extract the sequence
		coords=np.array(repeatAnnotDF[repeatAnnotDF['Name']==RNAid][['IndexStart','End_for_extraction']])
		sequence=np.array(list(repeat_genome_bases[coords[0][0]:coords[0][1]]))
		# Create data frame for storage 
		storageDF=pd.DataFrame()
		# Get RT stops
		RTpositions=np.loadtxt(filePath,delimiter='\t',dtype=int)
		# Histogram
		bins=range(coords[0][0],coords[0][1]+1,1) # Make sure bins are end coordinate inclusive
		hist,bins=np.histogram(RTpositions,bins=bins)
		width=0.7*(bins[1]-bins[0])
		center=(bins[:-1] + bins[1:])/2
		# Noramalize histogram frequencies by the total numer of RT stops 
		histPlot=np.array(hist,dtype=float)
		histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
		# Plot
		plt.subplot(plotDim,plotDim,i+1)
		plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
		plt.tick_params(axis='x',labelsize=2.5) 
		plt.tick_params(axis='y',labelsize=2.5)  
		plt.title('RT stops for %s: %s'%(RNAid,len(RTpositions)),fontsize=5)
		plt.xlim(coords[0][0],coords[0][1])
		# Extract sequence and record counts
		storageDF['Sequence']=sequence
		readsPerBase=np.array(list(hist))
		readsPerBaseNorm=np.array(list(histPlot))
		colName=sampleName+'RT_stops'
		storageDF[colName]=readsPerBase
		colNameNorm=sampleName+'RT_stops_norm'
		storageDF[colNameNorm]=readsPerBaseNorm
		# Layout
		plt.tight_layout()

def plot_rDNA(outfilepath,sampleName):
	repeatGenomeBuild=os.getcwd() + '/docs/repeat/repeatRNA.fa'
	repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
	repeat_genome_bases=repeat_genome[1]
	repeat_genome_size=len(repeat_genome[1])
	# Repeat index positions
	repeatAnnotation=os.getcwd() + '/docs/repeat/Hs_repeatIndex_positions.txt'
	repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
	repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
	# Python list extraction is not end index inclusive; to extract sequence, use end + 1.
	repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 

	ax=plt.subplot2grid((3,3),(0,0),colspan=3)
	RNAid='rDNA'
	repeatRNA_RTstops=glob.glob(outfilepath+'SourceData_RepeatRNA_*')
	filePath = [s for s in repeatRNA_RTstops if RNAid in s][0]
	coords=np.array(repeatAnnotDF[repeatAnnotDF['Name']==RNAid][['IndexStart','End_for_extraction']])
	rRNAstart=coords[0][0]
	sequence=np.array(list(repeat_genome_bases[coords[0][0]:coords[0][1]]))
	# Create data frame for storage 
	storageDF=pd.DataFrame()
	# Get RT stops
	RTpositions=np.loadtxt(filePath,delimiter='\t',dtype=int)
	# Histogram
	bins=range(coords[0][0],coords[0][1]+1,1) # Make sure bins are end coordinate inclusive
	hist,bins=np.histogram(RTpositions,bins=bins)
	width=0.7*(bins[1]-bins[0])
	center=(bins[:-1] + bins[1:])/2
	# Noramalize histogram frequencies by the total numer of RT stops 
	histPlot=np.array(hist,dtype=float)
	histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
	# Plot
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
	plt.tick_params(axis='x',labelsize=2.5) 
	plt.tick_params(axis='y',labelsize=2.5)  
	plt.title('RT stops for %s: %s'%(RNAid,len(RTpositions)),fontsize=5)
	plt.xlim(coords[0][0],coords[0][1])
	# Features of rDNA with respect to start of the bowtie index (index=0)
	start18s=3657
	end18s=5527
	start5s=6623
	end5s=6779
	start28s=7935
	end28s=12969
	# Overlay regions on plot 
	plt.axvspan(start18s+rRNAstart,end18s+rRNAstart,facecolor='g',alpha=0.5)
	plt.axvspan(start5s+rRNAstart,end5s+rRNAstart,facecolor='r',alpha=0.5)
	plt.axvspan(start28s+rRNAstart,end28s+rRNAstart,facecolor='b',alpha=0.5)

	# Generate histogram for transcribed region
	plt.subplot2grid((3,3),(1,0),colspan=3)
	# Generate histogram
	datarDNAOnly=RTpositions-rRNAstart
	bins=range((coords[0][0]-rRNAstart),(coords[0][1]-rRNAstart+1),1) # Make sure bins are end coordinate inclusive
	hist,bins=np.histogram(datarDNAOnly,bins=bins)
	width=0.7*(bins[1]-bins[0])
	center=(bins[:-1] + bins[1:])/2
	# Noramalize histogram frequencies by the total numer of RT stops 
	histPlot=np.array(hist,dtype=float)
	histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
	# Plot
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
	plt.tick_params(axis='x',labelsize=2.5) 
	plt.tick_params(axis='y',labelsize=2.5)  
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
	# Overlay regions
	plt.axvspan(start18s,end18s,facecolor='g',alpha=0.5)
	plt.axvspan(start5s,end5s,facecolor='r',alpha=0.5)
	plt.axvspan(start28s,end28s,facecolor='b',alpha=0.5)
	# Set x-axis to end position of rRNA locus
	rRNAend=13314
	plt.xlim(0,rRNAend)

	# Individual regions 
	plt.subplot2grid((3,3),(2,0),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='green',alpha=0.75)
	plt.xlim(start18s,end18s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('18s Region',fontsize=10)
	plt.subplot2grid((3,3),(2,1),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='red',alpha=0.75)
	plt.xlim(start5s,end5s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5) 
	plt.title('5.8s Region',fontsize=10)
	plt.subplot2grid((3,3),(2,2),colspan=1)
	plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.75)
	plt.xlim(start28s,end28s)
	plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
	plt.tick_params(axis='x',labelsize=5) 
	plt.tick_params(axis='y',labelsize=5)  
	plt.title('28s Region',fontsize=10)
	plt.tight_layout()

def plot_snoRNAdetail(outfilepath,sampleName):
	# Read through source files
	hcBOX_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_HCAbox.bed'
	SCA_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_SCARNA.bed'
	CDBOX_READS = outfilepath+'clipGenes_snoRNA_LowFDRreads_CDBOX.bed'
	readList=[SCA_READS,hcBOX_READS,CDBOX_READS]
	counts=[]
	for readFile in readList:
		counts.append(lineCount(readFile))
	
	totalCount=float(counts[0]+counts[1]+counts[2])

	ax=plt.subplot(2,2,1)
	# Check to make sure snoRNA reads recovered
	if totalCount > 0:
		inputVector=[float(counts[0])/totalCount,float(counts[1])/totalCount,float(counts[2])/totalCount]
		pie_wedges=ax.pie(inputVector,labels=["scaRNAs","hBox","CDbox"],labeldistance=1.1,autopct='%1.1f%%')
		plt.rcParams['font.size']=7.5
		for wedge in pie_wedges[0]:
			wedge.set_edgecolor('black')
			wedge.set_lw(1)
			plt.tight_layout()

	# Position files
	positionFiles=glob.glob(outfilepath+'*_genePosition_snoRNA')
	# Generate coverage histogram for each snoRNA type 
	snoNames=["C/D box","H/ACA box","scaRNA"]
	for j in range(0,len(positionFiles)):
		# Skip the first column, which has name of gene, and order is gene start / end and cluster start / end
		data=np.genfromtxt(positionFiles[j],usecols=(1,2,3,4,5),delimiter='\t',dtype='float')
		# Add exception handling for case in which no snoRNA reads come through
		try:
			# Ratio of read distance from start of gene to full gene length
			temp1=np.subtract(np.transpose(data[:,0]),np.transpose(data[:,2]))
			temp2=np.subtract(np.transpose(data[:,0]),np.transpose(data[:,1]))
			percentage=np.absolute(np.divide(temp1,temp2))
			# Correct any erronous values that are greater than one
			percentage[percentage > 1]=1
		except:
			print "No snoRNA reads of this type"
			percentage=np.zeros(1)
		# Plot histograms
		plt.subplot(2,2,j+2)
		bins=np.arange(0,1,0.01)
		hist,bins=np.histogram(percentage,bins=bins)
		hist=np.array(hist/float(percentage.size),dtype=float)
		width=0.7*(bins[1]-bins[0])
		center=(bins[:-1] + bins[1:])/2
		plt.bar(center, hist,align='center',width=width,color='blue',alpha=0.75)
		plt.tick_params(axis='x',labelsize=5) 
		plt.tick_params(axis='y',labelsize=5)  
		plt.xlabel('Fraction of gene body',fontsize=5)
		# For spacing the figures, only add ylabel to two bottom subplots
		if j+2 != 2:
			plt.ylabel('Normalized RT stop / bin',fontsize=5)
		plt.title(snoNames[j] + ' RT stops: %s'%str(percentage.size),fontsize=7.5)
		j += 1
		plt.tight_layout()

def plot_ncRNAs(outfilepath,sampleName):
	# Extract all files that list read start, stop as well as start,stop coordiantes for the gene.
	genePositionFiles=glob.glob(outfilepath+"*_genePosition")
	plotDim=math.ceil(math.sqrt(len(genePositionFiles)))
	j=1
	for i in range(0,len(genePositionFiles)):
		# Get the file name
		head, tail = os.path.split(genePositionFiles[i])
		geneTypeName = tail.split("_")[1]
		# Exclude certain gene names from the analysis; be sure to exclude snoRNAs since those are indepdently dealt with 
		if geneTypeName != 'IG' and geneTypeName != 'LRG' and geneTypeName != 'proteinCoding' and geneTypeName != 'TR' and geneTypeName != 'snRNA' and geneTypeName != 'rRNA':
			try:	
				# Skip the first column, which has name of gene, and order is gene start / end and cluster start / end
				data=np.loadtxt(genePositionFiles[i],dtype='float',usecols=range(1,5),delimiter='\t')
				clusterCenterPoint=np.zeros([data.shape[0],1])
				clusterCenterPoint=np.mean(data[:,[2,3]],axis=1)
				clusterCenterPoint=np.around(clusterCenterPoint,0)
				# Subtract cluster center from the gene start coordinate, and divide by gene length
				temp1=np.subtract(np.transpose(data[:,[0]]),clusterCenterPoint) # Must correct dimentions for row-wise subtraction 
				temp2=np.transpose(np.subtract(data[:,[0]],data[:,[1]]))
				percentage=np.absolute(np.divide(temp1,temp2))
				percentage[percentage > 1]=1 
			except:
				print "No data to read."
				percentage=np.zeros(1)	
			# Save source data for each ncRNA
			allReadCounts=outfilepath+'SourceData_ncRNAs_%s' %geneTypeName
			# Save array as a column vector
			np.savetxt(allReadCounts,percentage,delimiter='\t',fmt="%s")
			# Make subplot
			plt.subplot(3,4,j)
			bins=np.arange(0,1,0.01)
			hist,bins=np.histogram(percentage,bins=bins)
			# Noramlize histogram by total number of RT stops for the given ncRNA
			hist=np.array(hist/float(percentage.size),dtype=float)
			width=0.7*(bins[1]-bins[0])
			center=(bins[:-1] + bins[1:])/2
			plt.xlabel('Fraction of ncRNA gene body',fontsize=5)
			plt.ylabel('Normalized RT stop / bin',fontsize=5)
			plt.bar(center,hist,align='center',width=width,color='blue',alpha=0.75)
			plt.tick_params(axis='x',labelbottom='off') 
			plt.tick_params(axis='y',labelsize=5)  
			plt.title(geneTypeName+' RT stops: %s'%str(percentage.size),fontsize=5)
			plt.xlim(0,1)
			j += 1
			plt.tight_layout()

def runPlots(outfilepath,sampleName):

	fig1=plt.figure(1)
	plot_ReadAccounting(outfilepath,sampleName)
	plot_ReadsPerCluster(outfilepath,sampleName)
	plot_ClusterSizes(outfilepath,sampleName)
	plot_clusterBindingIntensity(outfilepath,sampleName)
	plot_readsBymRNAregion(outfilepath,sampleName)
	plot_BoundGeneTypes(outfilepath,sampleName)
	fig1.tight_layout()
	fig1.savefig(outfilepath+'Figure1.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig1.savefig(outfilepath+'Figure1.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

	fig2=plt.figure(2)
	plot_mRNAgeneBodyDist(outfilepath,sampleName)
	plot_geneBodyPartition(outfilepath,sampleName)
	fig2.tight_layout()
	fig2.savefig(outfilepath+'Figure2.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig2.savefig(outfilepath+'Figure2.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

	fig3=plt.figure(3)
	plot_repeatRNA(outfilepath,sampleName)
	fig3.tight_layout()
	fig3.savefig(outfilepath+'Figure3.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig3.savefig(outfilepath+'Figure3.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

	fig4=plt.figure(4)
	plot_rDNA(outfilepath,sampleName)
	fig4.tight_layout()
	fig4.savefig(outfilepath+'Figure4.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig4.savefig(outfilepath+'Figure4.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

	fig5=plt.figure(5)
	plot_snoRNAdetail(outfilepath,sampleName)
	fig5.tight_layout()
	fig5.savefig(outfilepath+'Figure5.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig5.savefig(outfilepath+'Figure5.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

	fig6=plt.figure(6)
	plot_ncRNAs(outfilepath,sampleName)
	fig6.tight_layout()
	fig6.savefig(outfilepath+'Figure6.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
	fig6.savefig(outfilepath+'Figure6.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

#############################################################################################################################

#############################################################################################################################

#############################################################################################################################

if __name__ == '__main__':

	# - iCLIP Parameters -

	# Chang data
	iCLIP5pBasesToTrim = 13 # Number of reads to trim from 5' end of clip reads.
	iCLIP3pBarcode = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' # Barcode for 3' end of CLIP reads.
	q=25 # Minimum quality score to keep during filtering.
	p=80 # Percentage of bases that must have quality > q during filtering.

	'''
	# Ule data
	iCLIP5pBasesToTrim = 9 # Number of reads to trim from 5' end of clip reads.
	iCLIP3pBarcode = 'AGATCGGAAGAGCGGTTCAG' # Barcode for 3' end of CLIP reads.
	phredQualitySystem='-Q33'
	'''

	# - Inputs - 
	global sampleName
	global outfilepath
	sampleName=sys.argv[1]
	infilepath = os.getcwd() + '/' + 'rawdata/' 
	outfilepath=os.getcwd() + '/results/' + sampleName + '/'

	# - Create log file - 
	global logFile
	logFile = outfilepath + "runLog"
	global logOpen
	logOpen = open(logFile, 'w')

	'''
	# - Get reads - 
	print "Processing sample %s" %sampleName
	unzippedreads = getFastq(infilepath,sampleName)

	logOpen.write("Processing sample: "+sampleName+'\n')
	# Create a new path to the output directory and move the files
	unzippedreads=moveFiles(unzippedreads,outfilepath)
	unzippedreads=changePath(unzippedreads,outfilepath)

	# - Trim 3p end of read files -
	print "Trim 3p"
	trimmedReads3p=trimReads3p(unzippedreads,iCLIP3pBarcode)

	# - Quality filter -
	print "Quality filter"
	filteredReads=qualityFilter(trimmedReads3p,q,p)

	# - Duplicate removal -
	print "Duplicate removal"
	nodupReads=dupRemoval(filteredReads)
	# Convert to fastq files
	nodupReadsFastq=fastaTofastq(nodupReads)

	# - Trim 5p end of read files -
	print "5p trim"
	trimmedReads5p=trimReads5p(nodupReadsFastq,iCLIP5pBasesToTrim)

	'''
	##########################
	#### START FOR MAPPING ###
	##########################
	# Trimmed reads
	# trimmedReads5p=glob.glob(outfilepath+'*_3ptrimmed_filter_nodupe_5ptrimmed.fastq')
	##########################
	
	'''

	# Internal dictionary for referencing genome index
	global genomeDict
	genomeDict = {}
	genomeDict['hg19']=1
	genomeDict['jfh1']=2
	genomeDict['h77']=3
	genomeDict['mm9']=4
	genomeDict['rep']=5

	# Run BOWTIE on repeat index and make bedFiles
	print "Run mapping to repeat index:"
	global genomeIndex
	genomeIndex='rep'
	mappedReads,unmappedReads=runBowtie(trimmedReads5p)
	mappedBedFiles=runSamtools(mappedReads)
	# Seperate reads by strand
	readsByStrand=seperateStrands(mappedBedFiles)
	# Isolate RT stops
	negativeRTstop=isolate5prime(modifyNegativeStrand(readsByStrand[0])) # negative strand reads
	positiveRTstop=isolate5prime(readsByStrand[1]) # positive strand reads 
	# Combine RT stops from replicates for each strand
	posMerged = outfilepath + sampleName + 'repeatRNA_positivereads.mergedRT'
	negMerged = outfilepath + sampleName + 'repeatRNA_negativereads.mergedRT'
	mergeRT(negativeRTstop,negMerged)
	mergeRT(positiveRTstop,posMerged)
	# Merge all
	negAndPosMerged = outfilepath + sampleName + 'repeatRNA_allreads.mergedRT'
	fileCat(negAndPosMerged,[posMerged,negMerged])

	# Run BOWTIE on hg19
	print "Run mapping to hg19:"
	genomeIndex='hg19'
	mappedReads,unmappedReads=runBowtie(unmappedReads)
	# Run Samtools
	mappedBedFiles=runSamtools(mappedReads)
	# Run Blacklist filter
	blacklistedBedFiles=runBlacklistRegions(mappedBedFiles)
	# Run repeat masker
	maskedBedFiles=runRepeatMask(blacklistedBedFiles)
	# Seperate reads by strand
	readsByStrand=seperateStrands(maskedBedFiles)
	# Isolate RT stops
	negativeRTstop=isolate5prime(modifyNegativeStrand(readsByStrand[0])) # negative strand reads
	positiveRTstop=isolate5prime(readsByStrand[1]) # positive strand reads 

	# Combine RT stops from replicates for each strand
	posMerged = outfilepath + sampleName + '_positivereads.mergedRT'
	negMerged = outfilepath + sampleName + '_negativereads.mergedRT'
	# Comment this out if starting from
	mergeRT(negativeRTstop,negMerged)
	mergeRT(positiveRTstop,posMerged)
	negAndPosMerged = outfilepath + sampleName + '_allreads.mergedRT.bed'
	fileCat(negAndPosMerged,[posMerged,negMerged])

	'''

	##########################
	#### START FOR CLIPPER, HOMER, AND PLOTTING ####
	##########################
	# Get sample name
	negAndPosMerged = outfilepath + sampleName + '_allreads.mergedRT.bed'
	##########################
	
	# --- Run CLIPPER ---
	print "Run CLIPPER:"
	CLIPPERio=runCLIPPER(negAndPosMerged)
	CLIPPERin=CLIPPERio[0]
	CLIPPERout=CLIPPERio[1]
	# Extract lowFDR reads based upon strand
	clipperStats=modCLIPPERout(CLIPPERin,CLIPPERout)
	CLIPPERlowFDR=clipperStats[0] # Low FDR reads returned filtred through CLIPper windows
	CLIPpeReadsPerCluster=clipperStats[1] # Number of reads per CLIPper cluster
	CLIPpeGeneList=clipperStats[2] # Gene names returned from the CLIPper file
	CLIPperOutBed=clipperStats[3] # CLIPper windows as a bed file
	# Sort CLIPPER clusters
	sortedClusters=sortCLIPClusters(CLIPperOutBed)
	# Make a bedGraph from all lowFDR CLIPper reads
	bedGraphCLIPout=makeBedGraph(CLIPPERlowFDR)
	# Get names of lowFDR genes of each gene type
	pathToGeneLists=getLowFDRGeneTypes(CLIPpeGeneList)
	# Extract all lowFDR reads of each gene type
	pathToReadLists=getLowFDRReadTypes(CLIPPERlowFDR,pathToGeneLists)

	# --- Filter and sort gene lists --- 
	print "Filter and sort gene lists:"
	proteinCodingReads = outfilepath+'clipGenes_proteinCoding_LowFDRreads.bed'
	proteinBedGraph=makeBedGraph(proteinCodingReads)
	proteinCentersBedGraph=makeBedGraph(getBedCenterPoints(proteinCodingReads))
	# Get cluster center points and filter snoRNAs and miRNAs from protein coding and snoRNAs
	filteredProteinCodingCenters=filterSnoRNAs(getBedCenterPoints(proteinCodingReads))
	filteredProteinCentersBedGraph=makeBedGraph(filteredProteinCodingCenters)
	lincRNAReads = outfilepath+'clipGenes_lincRNA_LowFDRreads.bed'
	filteredLincRNACenters=filterSnoRNAs(getBedCenterPoints(lincRNAReads))
	# Generate sorted gene lists
	pathToGeneLists=makeGeneLists(outfilepath,CLIPPERlowFDR)

	# --- Get intensity around cluster centers ---
	print "Get binding intensity around cluster centers"
	# Extract clipper input reads and cluster coordinates 
	clipperIN=outfilepath+sampleName+"_allreads.mergedRT_CLIPPERin.bed"
	CLIPclusters=outfilepath+sampleName+"_allreads.mergedRT_CLIP_clusters.bed"
	# Make bedgraph file for clipper reads 
	bedGraphCLIPin=makeBedGraph(clipperIN)
	# Extract the center coordinate of each clipper cluster 
	centerCoordinates=makeClusterCenter(CLIPclusters) 
	# Get histogram of read intensity around cluster center
	getClusterIntensity(bedGraphCLIPin,centerCoordinates)

	# --- Sort reads by gene body ---
	print "Intron and UTR analysis:"
	fivePreads,notFivePreads,CDSreads,notCDSreads,threePreads,notThreePreads=extractUTRs(filteredProteinCodingCenters)
	UTRgenelists=[fivePreads,CDSreads,threePreads]
	saveSourceData(UTRgenelists)

	# --- Detailed mRNA body analysis --- 
	print "Gene body analysis:"
	# Bed file with protein coding reads
	filteredProteinCoding = outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed'
	bedGraphProtein=makeBedGraph(filteredProteinCoding)
	# Generate gene-by-gene read coverage histogram 
	makeAvgGraph(bedGraphProtein)

	# --- Generate coverage histograms for ncRNAs ---
	print "Process non-coding RNAs:"
	readPositionFiles=[]
	for geneTypeReads in pathToGeneLists:
		# Do not perform this for protein coding or snoRNAs (performed using seperate gene reference file) 
		if 'snoRNA' not in geneTypeReads and 'proteinCoding' not in geneTypeReads:
			readPos=getReadDist(geneTypeReads)
			readPositionFiles=readPositionFiles+[readPos]
	partitionSnoRNAs(outfilepath)
	snoRNAfile = outfilepath + '/clipGenes_snoRNA_genePosition'
	fileList = [outfilepath+'/clipGenes_snoRNA_LowFDRreads_CDBOX_genePosition_snoRNA',outfilepath+'/clipGenes_snoRNA_LowFDRreads_HCAbox_genePosition_snoRNA',outfilepath+'/clipGenes_snoRNA_LowFDRreads_SCARNA_genePosition_snoRNA']
	fileCat(snoRNAfile,fileList)

	# --- Repeat RNA analysis ---
	print "Reapat analysis"
	negAndPosMerged = outfilepath + sampleName + 'repeatRNA_allreads.mergedRT'
	repeatMapped=parseRepeatMapped(negAndPosMerged)
	
	# ---- Plots ---

	# runPlots(outfilepath,sampleName)

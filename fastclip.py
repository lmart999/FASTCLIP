# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import cmath
import math
import sys
import numpy as np
import glob 
import subprocess
import re
from matplotlib_venn import venn2
import pandas as pd
from collections import defaultdict
from operator import itemgetter
import matplotlib as mpl
import matplotlib.pyplot as plt
import shutil
from optparse import OptionParser
mpl.rcParams['savefig.dpi'] = 2 * mpl.rcParams['savefig.dpi']

# <codecell>

global sampleName
global outfilepath
global logFile
global logOpen

### File name ###
sampleName=sys.argv[1]
infilepath=os.getcwd() + '/' + 'rawdata/'
outfilepath=os.getcwd() + '/results/%s/'%sampleName

# <codecell>

# Create log and start pipeline
logFile=outfilepath + "runLog"
logOpen=open(logFile, 'w')

# <codecell>

### Parameters ###
iCLIP3pBarcode='AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG' # Barcode sequence to trim from reads.
q=25 # Minimum quality score to keep during filtering.
p=80 # Percentage of bases that must have quality > q during filtering.    
iCLIP5pBasesToTrim=13 # Number of reads to trim from 5' end of clip reads.
k='1' # k=N distinct, valid alignments for each read in bt2 mapping.
threshold=7 # Sum of RT stops (for both replicates) required to keep file. 
expand=15 # Bases to expand around RT position after RT stops are merged.
repeat_index=os.getcwd() + '/docs/repeat/rep' # bt2 index for repeat RNA.
repeatGenomeBuild=os.getcwd()+'/docs/repeat/repeatRNA.fa' # Sequence of repeat index.
repeatAnnotation=os.getcwd()+'/docs/repeat/Hs_repeatIndex_positions.txt' # Repeat annotation file.
start18s=3657
end18s=5527
start5s=6623
end5s=6779
start28s=7935
end28s=12969
rRNAend=13314
threshold_rep=0 # RT stop threshold for repeat index.
index=os.getcwd() + '/docs/hg19/hg19' # bt2 index for mapping.
index_tag='hg19' # Name of bt2 index.
genomeFile=os.getcwd()+'/docs/human.hg19.genome' # Genome file for bedGraph, etc.
genomeForCLIPper='-shg19' # Parameter for CLIPper.
blacklistregions=os.getcwd()+'/docs/wgEncodeDukeMapabilityRegionsExcludable.bed' # Blacklist masker.
repeatregions=os.getcwd()+'/docs/repeat_masker.bed' # Repeat masker.
geneAnnot=glob.glob(os.getcwd()+'/docs/genes_types/*') # List of genes by type.
snoRNAmasker=os.getcwd()+'/docs/snoRNA_reference/snoRNAmasker_formatted_5pExtend.bed' # snoRNA masker file.
miRNAmasker=os.getcwd()+'/docs/miR_sort_clean.bed' # miRNA masker file.
fivePUTRBed=os.getcwd()+'/docs/5pUTRs_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
threePUTRBed=os.getcwd()+'/docs/3pUTRs_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
cdsBed=os.getcwd()+'/docs/Exons_Ensbl_sort_clean_uniq.bed' # UTR annotation file.
utrFile=os.getcwd()+'/docs/hg19_ensembl_UTR_annotation.txt' # UTR annotation file.
genesFile=os.getcwd()+'/docs/hg19_ensembl_genes.txt' # Gene annotation file.
sizesFile=os.getcwd()+'/docs/hg19.sizes' # Genome sizes file. 
snoRNAindex=os.getcwd()+'/docs/snoRNA_reference/sno_coordinates_hg19_formatted.bed' # snoRNA coordinate file.
CLIPPERoutNameDelim='_' # Delimiter that for splitting gene name in the CLIPper windows file.

# <codecell>

import datetime
now=datetime.datetime.now()
logOpen.write("Timestamp:%s\n"%str(now))
logOpen.write("\n###Parameters used###\n")
logOpen.write("3' barcode:%s\n'"%iCLIP3pBarcode)
logOpen.write("Minimum quality score (q):%s\n"%q)
logOpen.write("Percentage of bases with > q:%s\n"%p)
logOpen.write("5' bases to trim:%s\n'"%iCLIP5pBasesToTrim)
logOpen.write("k distinct, valid alignments for each read in bt2 mapping:%s\n"%k)
logOpen.write("Threshold for minimum number of RT stops:%s\n"%threshold)
logOpen.write("Bases for expansion around conserved RT stops:%s\n"%expand)
logOpen.write("\n\n\n")

# <codecell>

print "Processing sample %s" %(sampleName)
logOpen.write("Processing sample: "+sampleName+'\n')
read1=infilepath+sampleName+'_R1.fastq'
read2=infilepath+sampleName+'_R2.fastq'
unzippedreads=[read1,read2]

# <codecell>

def trimReads3p(unzippedreads,adapter3p):
    # Usage: Trims a specified adapter sequence from the 3p end of the reads.
    # Input: List of fastq files.
    # Output: List of 3p trimmed files.
    trimparam='-a'+adapter3p # Adapter string
    trimmedReads=[]
    try:
        for inread in unzippedreads:
            outread=inread.replace("rawdata/", "results/%s/"%sampleName)
            outread=outread.replace(".fastq", "_3ptrimmed.fastq")
            process=subprocess.Popen(['fastx_clipper',trimparam,'-n','-l33','-Q33','-i',inread,'-o',outread],stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            logOpen.write("Trim 3p end of reads.\n")
            logOpen.write("Stdout: %s.\n"%stdout)
            logOpen.write("Stderr: %s.\n"%stderr)
            trimmedReads=trimmedReads+[outread]
        return trimmedReads
    except:
        logOpen.write("Problem with 3p trimming.\n")
        print "Problem with 3p trimming."

print "Trim 3p adapter from reads."
trimmedReads3p=trimReads3p(unzippedreads,iCLIP3pBarcode)

# <codecell>

def qualityFilter(trim3pReads,q,p):
    # Usage: Filters reads based upon quality score.
    # Input: List of fastq file names as well as the quality paramters p and q.
    # Output: List of modified fastq file names.
    qualityparam='-q'+str(q)
    percentrageparam='-p'+str(p)
    filteredReads=[]
    try:
        for inread in trim3pReads:
            outread=inread.replace(".fastq", "_filter.fastq")
            process=subprocess.Popen(['fastq_quality_filter',qualityparam,percentrageparam,'-Q33','-i',inread,'-o',outread],stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
            stdout, stderr=process.communicate()
            logOpen.write("Perform quality filtering.\n")
            logOpen.write("Stdout: %s.\n"%stdout)
            logOpen.write("Stderr: %s.\n"%stderr)
            filteredReads=filteredReads+[outread]
        return filteredReads
    except:
        logOpen.write("Problem with quality filter.\n")
        print "Problem with quality filter."

print "Perform quality filtering."
filteredReads=qualityFilter(trimmedReads3p,q,p)

# <codecell>

def dupRemoval(filteredReads):
    # Usage: Removes duplicate reads.
    # Input: List of fastq file names.
    # Output: List of reads in FASTA format.
    program=os.getcwd() + '/bin/fasta_to_fastq.pl'
    noDupes=[]
    try:
        for inread in filteredReads:
            outread=inread.replace(".fastq","_nodupe.fasta")
            process=subprocess.Popen(['fastx_collapser','-Q33','-i',inread,'-o',outread],stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
            stdout, stderr=process.communicate()
            logOpen.write("Perform duplicate removal.\n")
            logOpen.write("Stdout: %s.\n"%stdout)
            logOpen.write("Stderr: %s.\n"%stderr)
            fastqOut=outread.replace('.fasta', '.fastq') # fastx_collapser returns fasta files, which are then converted to fastq.
            outfh=open(fastqOut, 'w')
            process=subprocess.Popen(['perl',program,outread],stdout=outfh)
            process.communicate() # Wait for the process to complete.
            os.remove(outread) # Remove the remaining .fasta file.
            noDupes=noDupes+[fastqOut]
        return noDupes
    except:
        logOpen.write("Problem with duplicate removal.\n")
        print "Problem with duplicate removal."
        
print "Perform duplicate removal."
nodupReads=dupRemoval(filteredReads)

# <codecell>

def trimReads5p(nodupes,n):
    # Usage: Trims a specified number of bases from the 5' end of each read.
    # Input: List of fastq files.
    # Output: List of 5p trimmed files.
    trimparam='-f'+str(n)
    trimmedReads=[]
    try:
        for inread in nodupes:
            outread=inread.replace(".fastq", "_5ptrimmed.fastq")
            process=subprocess.Popen(['fastx_trimmer', trimparam, '-Q33', '-i', inread,'-o',outread],stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
            stdout, stderr=process.communicate()
            logOpen.write("Perform 5' barcode trimming.\n")
            logOpen.write("Stdout: %s.\n"%stdout)
            logOpen.write("Stderr: %s.\n"%stderr)
            trimmedReads=trimmedReads+[outread]
        return trimmedReads
    except:
        logOpen.write("Problem with 5' barcode trimming.\n")
        print "Problem with 5' barcode trimming."

print "Perform 5' barcode trimming."
trimmedReads5p=trimReads5p(nodupReads,iCLIP5pBasesToTrim)

# <codecell>

def runBowtie(fastqFiles,index,index_tag):
    # Usage: Read mapping to reference.
    # Input: Fastq files of replicate trimmed read files.
    # Output: Path to samfile for each read.
    program='bowtie2'
    mappedReads=[]
    unMappedReads=[]
    try:
        for infastq in fastqFiles:
            outfile=infastq.replace(".fastq","_mappedTo%s.sam"%index_tag)
            unmapped=infastq.replace(".fastq","_notMappedTo%s.fastq"%index_tag)
            process=subprocess.Popen([program,'-x',index,'-k',k,'-U',infastq,'--un',unmapped,'-S',outfile],stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
            stdout,stderr=process.communicate()
            logOpen.write("Perform mapping to %s index.\n"%index_tag)
            logOpen.write("Stdout: %s.\n"%stdout)
            logOpen.write("Stderr: %s.\n"%stderr)        
            mappedReads = mappedReads + [outfile]
            unMappedReads = unMappedReads + [unmapped]
        return (mappedReads,unMappedReads)
    except:
        logOpen.write("Problem with mapping.\n")
        print "Problem with mapping."

print "Run mapping to repeat index."  
mappedReads_rep,unmappedReads_rep=runBowtie(trimmedReads5p,repeat_index,'repeat')

# <codecell>

def runSamtools(samfiles):
    # Usage: Samfile processing.
    # Input: Sam files from Bowtie mapping.
    # Output: Sorted bedFiles.
    program = 'samtools'
    program2 = 'bamToBed'
    outBedFiles=[]
    try:
        for samfile in samfiles:
            bamfile = samfile.replace('.sam','.bam')  
            proc = subprocess.Popen( [program,'view','-bS','-o', bamfile, samfile]) 
            proc.communicate()
            bamfile_sort = bamfile.replace('.bam','_sorted') 
            proc2 = subprocess.Popen([program,'sort',bamfile, bamfile_sort])
            proc2.communicate()
            bedFile = bamfile_sort.replace('_sorted', '_withDupes.bed') 
            outfh = open(bedFile,'w')
            proc3 = subprocess.Popen( [program2,'-i', bamfile_sort+'.bam'],stdout=outfh)
            proc3.communicate()
            outBedFiles=outBedFiles+[bedFile]   
        return outBedFiles
    except:
        logOpen.write("Problem with samtools.\n")
        print "Problem with samtools."

print "Run samtools."
logOpen.write("Run samtools.\n")
mappedBedFiles_rep=runSamtools(mappedReads_rep)

# <codecell>

def seperateStrands(mappedReads):
	# Usage: Seperate positive and negative strands.
	# Input: Paths to two bed files from Samtools.
	# Output: Paths to bed files isolated by strand.
    negativeStrand=[]
    positiveStrand=[]
    for mapFile in mappedReads:
        with open(mapFile, 'r') as infile:
            neg_strand=mapFile.replace('.bed','_neg.bed')
            pos_strand=mapFile.replace('.bed','_pos.bed')	
            neg = open(neg_strand, 'w')
            pos = open(pos_strand, 'w')
            negativeStrand=negativeStrand+[neg_strand]
            positiveStrand=positiveStrand+[pos_strand]
            for line in infile:	
                if str(line.strip().split('\t')[5]) == '-':
                    neg.write(line)
                elif str(line.strip().split('\t')[5]) == '+':
                    pos.write(line)
    return (negativeStrand,positiveStrand)

def modifyNegativeStrand(negativeStrandReads):
    # Usage: For negative stranded reads, ensure 5' position (RT stop) is listed first.
    # Input: Bed file paths to all negative stranded.
    # Output: Paths to modified bed files.
    negativeStrandEdit=[]
    for negativeRead in negativeStrandReads:
        neg_strand_edited=negativeRead.replace('_neg.bed','_negEdit.bed')
        negativeStrandEdit=negativeStrandEdit+[neg_strand_edited]
        neg_edit = open(neg_strand_edited, 'w')
        with open(negativeRead, 'r') as infile:
            for line in infile:	
                chrom,start,end,name,quality,strand=line.strip().split('\t')
                neg_edit.write('\t'.join((chrom,end,str(int(end)+30),name,quality,strand))+'\n')
    return negativeStrandEdit

def isolate5prime(strandedReads):
	# Usage: Isolate only the Chr, 5' position (RT stop), and strand.
	# Input: Bed file paths to strand seperated reads.
	# Output: Paths RT stop files.
    RTstops=[]
    for reads in strandedReads:
        RTstop=reads.replace('.bed','_RTstop.bed')
        f = open(RTstop, 'w')
        with open(reads, 'r') as infile:
            RTstops=RTstops+[RTstop]
            for line in infile:	
                chrom,start,end,name,quality,strand=line.strip().split('\t')
                f.write('\t'.join((chrom,start,strand))+'\n')
    return RTstops

print "RT stop isolation (repeat)."
logOpen.write("RT stop isolation (repeat).\n")
readsByStrand_rep=seperateStrands(mappedBedFiles_rep)
negativeRTstop_rep=isolate5prime(modifyNegativeStrand(readsByStrand_rep[0])) 
positiveRTstop_rep=isolate5prime(readsByStrand_rep[1]) 

# <codecell>

def fileCat(destinationFile,fileList):
    f = open(destinationFile, "w")
    for tempfile in fileList:
        readfile = open(tempfile, "r")
        f.write(readfile.read())
        readfile.close()
    f.close()

def RTcounts(RTfile):
    posRT_R1=pd.DataFrame(pd.read_table(RTfile,index_col=None,header=None,sep='\t'))
    posRT_R1.columns=['Chr','Start','Strand']
    cts=posRT_R1.groupby(['Chr','Start']).size()
    return cts

def mergeRT(RTstopFiles,outfilename,threshold,expand,strand):
    # Usage: Merge RT stops between replicates and keep only those positions that exceed threshold.
    # Input: Files with RT stops for each replicate, outfile, threshold, strand, and bases to expand around RT stop.
    # Output: None. Writes merged RT stop file.
    cts_R1=RTcounts(RTstopFiles[0])
    cts_R2=RTcounts(RTstopFiles[1])
    m=pd.concat([cts_R1,cts_R2],axis=1,join='inner')
    m.columns=['Rep_1','Rep_2']
    m['Sum']=m['Rep_1']+m['Rep_2']
    m_filter=m[m['Sum']>threshold]
    f = open(outfilename, 'w')
    for i in m_filter.index:
        chrom=i[0]
        RT=i[1]
        count=m_filter.loc[i,'Sum']
        if RT > expand:
            read='\t'.join((chrom,str(int(RT)-expand),str(int(RT)+expand),'CLIPread','255',strand))+'\n'
        else:
            read='\t'.join((chrom,str(int(RT)),str(int(RT)+expand),'CLIPread','255',strand))+'\n'
        f.write(read*(count))

print "Merge RT stops."
logOpen.write("Merge RT stops.\n")
posMerged=outfilepath+sampleName+'_repeat_positivereads.mergedRT'
strand='+'
mergeRT(positiveRTstop_rep,posMerged,threshold_rep,expand,strand)
negMerged=outfilepath+sampleName+'_repeat_negativereads.mergedRT'
strand='-'
mergeRT(negativeRTstop_rep,negMerged,threshold_rep,expand,strand)
negAndPosMerged=outfilepath+sampleName+'_threshold=%s'%threshold_rep+'_repeat_allreads.mergedRT.bed'
fileCat(negAndPosMerged,[posMerged,negMerged])

# <codecell>

print "Run mapping to %s."%index_tag
mappedReads,unmappedReads=runBowtie(unmappedReads_rep,index,index_tag)

# <codecell>

print "Run samtools."
logOpen.write("Run samtools.\n")
mappedBedFiles=runSamtools(mappedReads)

# <codecell>

def runRepeatMask(mappedReads,repeatregions):
	# Usage: Remove repeat regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools) and blastlist regions removed.
	# Output: Bedfile with repeat regions removed.
	program='intersectBed'
	masked=[]
	try:
		for bedIn in mappedReads:
			noRepeat=bedIn.replace('.bed','_noRepeat.bed')
			outfh=open(noRepeat, 'w')
			proc=subprocess.Popen([program,'-a',bedIn,'-b',repeatregions,'-v','-s'],stdout=outfh)
			proc.communicate()
			outfh.close()
			masked=masked+[noRepeat]
		return (masked)
	except:
		print "Problem with repeat masking."
        logOpen.write("Problem with repeat masking.\n")
        
def runBlacklistRegions(mappedReads,blacklistregions):
	# Usage: Remove blacklisted regions from bedfile following mapping.
	# Input: .bed file after mapping (duplicates removed by samtools).
	# Output: Bedfile with blacklisted regions removed.
	program='intersectBed'
	blackListed=[]
	try:
		for bedIn in mappedReads:
			noBlacklist=bedIn.replace('.bed','_noBlacklist.bed')
			outfh=open(noBlacklist, 'w')
			proc=subprocess.Popen([program,'-a',bedIn,'-b',blacklistregions,'-v'],stdout=outfh)
			proc.communicate()
			outfh.close()
			blackListed=blackListed+[noBlacklist]
		return (blackListed)
	except:
		print "Problem with blacklist."
        logOpen.write("Problem with blacklist.\n")
        
print "Run repeat and blacklist region masker."
logOpen.write("Run repeat and blacklist masker.\n")
blacklistedBedFiles=runBlacklistRegions(mappedBedFiles,blacklistregions)
maskedBedFiles=runRepeatMask(blacklistedBedFiles,repeatregions)

# <codecell>

print "RT stop isolation."
logOpen.write("RT stop isolation.\n")
readsByStrand=seperateStrands(maskedBedFiles)
negativeRTstop=isolate5prime(modifyNegativeStrand(readsByStrand[0])) 
positiveRTstop=isolate5prime(readsByStrand[1]) 

print "Merge RT stops."
logOpen.write("Merge RT stops.\n")
posMerged=outfilepath+sampleName+'_%s_positivereads.mergedRT'%index_tag
strand='+'
mergeRT(positiveRTstop,posMerged,threshold,expand,strand)
negMerged=outfilepath+sampleName+'_%s_negativereads.mergedRT'%index_tag
strand='-'
mergeRT(negativeRTstop,negMerged,threshold,expand,strand)
negAndPosMerged=outfilepath+sampleName+'_threshold=%s'%threshold+'_%s_allreads.mergedRT.bed'%index_tag
fileCat(negAndPosMerged,[posMerged,negMerged])

# <codecell>

def runCLIPPER(RTclusterfile,genome,genomeFile):
    # Useage: Process the mergedRT file and pass through CLIPper FDR script.
    # Input: Merged RT file.
    # Output: CLIPper input (.bed) file and output file.
    program='bedToBam'
    program2='samtools'
    program3='bamToBed'
    program4='clipper'
    
    bamfile=RTclusterfile.replace('.bed','.bam')  
    outfh=open(bamfile, 'w')
    proc=subprocess.Popen([program,'-i',RTclusterfile,'-g',genomeFile],stdout=outfh)
    proc.communicate()
    
    bamfile_sort=bamfile.replace('.bam','.srt')
    proc2=subprocess.Popen([program2,'sort',bamfile,bamfile_sort])
    proc2.communicate()
    
    bamfile_sorted=bamfile_sort+'.bam'
    mapStats=bamfile_sorted.replace('.srt.bam','.mapStats.txt') 
    outfh=open(mapStats, 'w')
    proc3=subprocess.Popen([program2,'flagstat',bamfile_sorted],stdout=outfh)
    proc3.communicate()
    
    proc4=subprocess.Popen([program2,'index',bamfile_sorted])
    proc4.communicate()
    
    CLIPPERin=bamfile_sorted.replace('.srt.bam','_CLIPPERin.bed') 
    outfh=open(CLIPPERin, 'w')
    proc5=subprocess.Popen([program3,'-i',bamfile_sorted],stdout=outfh)
    proc5.communicate()
    
    CLIPPERout=CLIPPERin.replace('_CLIPPERin.bed','_CLIP_clusters') 
    proc6=subprocess.Popen([program4,'--bam',bamfile_sorted,genome,'--outfile=%s'%CLIPPERout],)
    proc6.communicate()
    outfh.close()
    
    return (CLIPPERin,CLIPPERout)

def makeGeneNameDict(fi):
    # Usage: Make a dictionary that maps RT stop to gene name.
    # Input: File path to intersected CLIPper windows and input RT stop coordinates.
    # Output Dictionary mapping RT stop to name.
    nameDict={}
    with open(fi, 'r') as infile:
        for read in infile:
            elementList=read.strip().split('\t')
            RT_id='_'.join((elementList[0],elementList[1],elementList[2],elementList[5]))
            if RT_id not in nameDict:
                geneName=elementList[9].strip().split(CLIPPERoutNameDelim)[0]
                nameDict[RT_id]=geneName
    return nameDict

def modCLIPPERout(CLIPPERin,CLIPPERout):
    # Usage: Process the CLIPper output and isolate lowFDR reads based upon CLIPper windows.
    # Input: .bed file passed into CLIPper and the CLIPper windows file.
    # Output: Low FDR reads recovered using the CLIPer windows file, genes per cluster, gene list of CLIPper clusters, and CLIPper windows as .bed.
    program='intersectBed'
    CLIPperOutBed=CLIPPERout+'.bed'
    CLIPpeReadsPerCluster=CLIPPERout+'.readsPerCluster'
    CLIPpeGeneList=CLIPPERout+'.geneNames'
    f = open(CLIPperOutBed,'w')
    g = open(CLIPpeReadsPerCluster,'w')
    h = open(CLIPpeGeneList,'w')
    with open(CLIPPERout,'r') as infile:
        for line in infile:	
            try:
                # Note that different versions on CLIPper will report the gene name differently. So, we must handle this.
                chrom,start,end,name,stats,strand,start_2,end_2 = line.strip().split('\t')
                if CLIPPERoutNameDelim=='_':
                    readPerCluster=name.strip().split(CLIPPERoutNameDelim)[2]
                else: 
                    readPerCluster=(name.strip().split(CLIPPERoutNameDelim)[1]).split('_')[2]
                geneName=name.strip().split(CLIPPERoutNameDelim)[0]
                f.write('\t'.join((chrom,start,end,name,stats,strand))+'\n')
                g.write((readPerCluster+'\n'))
                h.write((geneName+'\n'))
            except:
                print ""
    f.close()
    g.close()
    h.close()
    
    # Intersect input reads with the CLIPper windows, report full result for both, include strand, do not duplicate reads from -a if they interset with multiple windows.
    clusterWindowInt=CLIPperOutBed.replace('.bed','_fullClusterWindow.bed')
    outfh=open(clusterWindowInt,'w')
    proc=subprocess.Popen([program,'-a',CLIPPERin,'-b',CLIPperOutBed,'-wa','-wb','-s'],stdout=outfh)
    proc.communicate()
    outfh.close()
    
    # Use the full window intersection to make a dictionary mapping RTstop to gene name.
    nameDict=makeGeneNameDict(clusterWindowInt)
    
    # Intersect input reads with CLIPper windows, but only report one intersection per read (as reads can overlap with multiple windows).
    clusterWindowIntUniq=CLIPperOutBed.replace('.bed','_oneIntPerRead.bed')
    outfh=open(clusterWindowIntUniq,'w')
    proc=subprocess.Popen([program,'-a',CLIPPERin,'-b',CLIPperOutBed,'-wa','-s','-u'],stdout=outfh)
    proc.communicate()
    outfh.close()
    
    # Process the uniquly intersected RT stops by adding gene name.
    CLIPPERlowFDR=CLIPperOutBed.replace('.bed','_lowFDRreads.bed')
    outfh=open(CLIPPERlowFDR,'w')
    with open(clusterWindowIntUniq, 'r') as infile:
        for read in infile:
            bed=read.strip().split('\t')
            RT_id='_'.join((bed[0],bed[1],bed[2],bed[5]))
            geneName=nameDict[RT_id]
            outfh.write('\t'.join((bed[0],bed[1],bed[2],geneName,bed[4],bed[5],'\n')))
    outfh.close()
    infile.close()
    
    return (CLIPPERlowFDR,CLIPpeReadsPerCluster,CLIPpeGeneList,CLIPperOutBed)

print "Run CLIPper."
logOpen.write("Run CLIPper.\n")
CLIPPERio=runCLIPPER(negAndPosMerged,genomeForCLIPper,genomeFile)
CLIPPERin=CLIPPERio[0]
CLIPPERout=CLIPPERio[1]
clipperStats=modCLIPPERout(CLIPPERin,CLIPPERout)
CLIPPERlowFDR=clipperStats[0] # Low FDR reads returned filtred through CLIPper windows
CLIPpeReadsPerCluster=clipperStats[1] # Number of reads per CLIPper cluster
CLIPpeGeneList=clipperStats[2] # Gene names returned from the CLIPper file
CLIPperOutBed=clipperStats[3] # CLIPper windows as a bed file

# <codecell>

def getBedCenterPoints(inBed):
    # Usage: Obtain ceter coordiantes of bedFile.
    # Input: BedFile.
    # Output: Center coodinates returned.
    outBed=inBed.replace('.bed','_centerCoord.bed')	
    f=open(outBed, 'w')
    with open(inBed, 'r') as infile:
        for line in infile:	
            elementList=line.strip().split('\t')
            f.write('\t'.join((elementList[0],str(int(elementList[1])+expand),str(int(elementList[1])+expand+1),elementList[3],elementList[4],elementList[5],'\n')))
    f.close()
    return outBed

def cleanBedFile(inBed):
    # Usage: Sort and recover only first 6 fields from a bed file.
    # Input: BedFile.
    # Output: Sorted bedFile with correct number of fields.
    program='sortBed'
    CLIPperOutBed=inBed.replace('.bed','_cleaned.bed')	
    sortedBed=CLIPperOutBed.replace('_cleaned.bed','_cleaned_sorted.bed')
    f=open(CLIPperOutBed, 'w')
    with open(inBed, 'r') as infile:
        for line in infile:	
            elementList=line.strip().split('\t')
            f.write('\t'.join((elementList[0],elementList[1],elementList[2],elementList[3],elementList[4],elementList[5],'\n')))
    f.close()
    outfh=open(sortedBed, 'w')
    proc=subprocess.Popen([program, '-i', CLIPperOutBed],stdout=outfh)
    proc.communicate()
    outfh.close()
    return sortedBed

def makeBedGraph(lowFDRreads,sizesFile):
    # Usage: From a bedFile, generate a bedGraph and bigWig.
    # Input: BedFile.
    # Output: BedGraph file.
    program='genomeCoverageBed'
    program2=os.getcwd() + '/bin/bedGraphToBigWig'
    cleanBed=cleanBedFile(lowFDRreads)
    outname=cleanBed.replace('.bed','.bedgraph')
    outname2=cleanBed.replace('.bed','.bw')
    outfh=open(outname,'w')
    proc=subprocess.Popen([program,'-bg','-split','-i',cleanBed,'-g',sizesFile],stdout=outfh)
    proc.communicate()
    outfh2=open(outname2,'w')
    proc2=subprocess.Popen([program2,outname,sizesFile,outname2],stdout=subprocess.PIPE)
    proc2.communicate()
    return outname

print "Make bedGraph"
logOpen.write("Make bedGraph.\n")
bedGraphCLIPout=makeBedGraph(CLIPPERlowFDR,genomeFile)
CLIPPERlowFDRcenters=getBedCenterPoints(CLIPPERlowFDR)
allLowFDRCentersBedGraph=makeBedGraph(CLIPPERlowFDRcenters,genomeFile)

# <codecell>

def filterSnoRNAs(proteinCodingReads,snoRNAmasker,miRNAmasker):
    # Usage: Filter snoRNA and miRNAs from protein coding reads.
    # Input: .bed file with protein coding reads.
    # Output: snoRNA and miR filtered .bed file.
    program='intersectBed'
    proteinWithoutsnoRNAs=proteinCodingReads.replace('.bed','_snoRNAremoved.bed')
    proteinWithoutmiRNAs=proteinWithoutsnoRNAs.replace('.bed','_miRNAremoved.bed')
    outfh=open(proteinWithoutsnoRNAs, 'w')
    proc=subprocess.Popen([program,'-a',proteinCodingReads,'-b',snoRNAmasker,'-v','-s'],stdout=outfh)
    proc.communicate()
    outfh.close()
    outfh=open(proteinWithoutmiRNAs, 'w')
    proc=subprocess.Popen([program,'-a',proteinWithoutsnoRNAs,'-b',miRNAmasker,'-v','-s'],stdout=outfh)
    proc.communicate()
    outfh.close()
    return (proteinWithoutmiRNAs)

def getLowFDRReadTypes(CLIPPERlowFDR,pathToGeneLists):
    # Usage: Given a list of genes, return all reads for the associated genes.
    # Input: Gene list and the path to lowFDR read file.
    # Output: List of reads assocaited with the given genes.
    lowFDRgenelist=[]
    for path in pathToGeneLists:
        outfile=path+'_LowFDRreads.bed'
        proc=subprocess.Popen('grep -F -f %s %s > %s'%(path,CLIPPERlowFDR,outfile),shell=True)
        proc.communicate()
        return_code=proc.wait() # *** Remove later. ***
        lowFDRgenelist=lowFDRgenelist+[outfile]
    return lowFDRgenelist

def compareLists(list1,list2,outname):
    # Usage: Compare gene lists and output matches to the file. 
    # Input: Two gene lists.
    # Output: Path file containing the matching genes.
    f=open(list1,'r')
    g=open(list2,'r')
    commonGenes=set(f.readlines()) & set(g.readlines())
    geneCategory=outname.split('.')[1]
    outputName=outfilepath+'clipGenes_'+geneCategory
    outfh=open(outputName,'w')
    for gene in commonGenes:
        outfh.write(gene)
    outfh.close()
    return outputName

def getLowFDRGeneTypes(CLIPpeGeneList,geneAnnot):
    # Usage: Get all genes listed under each type, compare to CLIPper targets.
    # Input: .bed file passed into CLIPper and the CLIPper windows file.
    # Output: Path to file containing all CLIPper genes of each type.
    geneTypes=[]
    for genepath in geneAnnot:
        lowFDRgenes=compareLists(CLIPpeGeneList,genepath,os.path.split(genepath)[1])
        geneTypes=geneTypes+[lowFDRgenes]
    return geneTypes

print "Partition reads by type."
logOpen.write("Partition reads by type.\n")

pathToGeneLists=getLowFDRGeneTypes(CLIPpeGeneList,geneAnnot)
pathToReadLists=getLowFDRReadTypes(CLIPPERlowFDR,pathToGeneLists)

proteinCodingReads=outfilepath+'clipGenes_proteinCoding_LowFDRreads.bed'
proteinBedGraph=makeBedGraph(proteinCodingReads,genomeFile)
filteredProteinCodingCenters=filterSnoRNAs(getBedCenterPoints(proteinCodingReads),snoRNAmasker,miRNAmasker)
filteredProteinCentersBedGraph=makeBedGraph(filteredProteinCodingCenters,genomeFile)

lincRNAReads=outfilepath+'clipGenes_lincRNA_LowFDRreads.bed'
filteredLincRNACenters=filterSnoRNAs(getBedCenterPoints(lincRNAReads),snoRNAmasker,miRNAmasker)

# <codecell>

# --- # 

# <codecell>

def sortFilteredBed(bedFile):
    bf=pd.DataFrame(pd.read_table(bedFile,header=None))
    bf.columns=['Chr','Start','Stop','CLIPper_name','Q','Strand']
    geneCounts=countHitsPerGene(bf)
    return geneCounts

def countHitsPerGene(bf):
    # *** THIS MAY DEPEND UPON THE VERSION OF CLIPPER USED ***
    bf['geneName']=bf['CLIPper_name'].apply(lambda x: x.split('_')[0])
    geneCounts=bf.groupby('geneName').size()
    geneCounts.sort(ascending=False)
    return geneCounts

def getSnoRNAreads(CLIPPERlowFDRcenters,snoRNAindex):
    program='intersectBed'		
    bedFile=outfilepath+'clipGenes_snoRNA_LowFDRreads.bed'
    outfh=open(bedFile, 'w')
    proc=subprocess.Popen([program,'-a',CLIPPERlowFDRcenters,'-b',snoRNAindex,'-s','-wa','-wb'],stdout=outfh)
    proc.communicate()
    outfh.close()	
    return bedFile

def countSnoRNAs(bedFile_sno):
    bf=pd.DataFrame(pd.read_table(bedFile_sno,header=None))
    bf.columns=['Chr','Start','End','CLIPper_name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
    geneCounts=bf.groupby('name_snoRNA').size()
    geneCounts.sort(ascending=False)
    return geneCounts

def countRemainingGeneTypes(remaining):
    for bedFile in remaining:
        try:
            bf=pd.DataFrame(pd.read_table(bedFile,header=None))
            bf.columns=['Chr','Start','End','ReadName','Q','Strand','CLIPper_winChr','CLIPper_winStart','CLIPper_winEmd','CLIPper_winaName','CLIPper_winP','CLIPper_winStrand']
            # *** THIS MAY DEPEND UPON THE VERSION OF CLIPPER USED ***
            bf['geneName']=bf['CLIPper_winaName'].apply(lambda x: x.split('_')[0])
            geneCounts=bf.groupby('geneName').size()
            geneCounts.sort(ascending=False) 
                        
            head,fname=os.path.split(bedFile)
            geneType=fname.split("_")[1]
            outfilepathToSave=outfilepath+'/PlotData_ReadsPerGene_%s'%geneType
            geneCounts.to_csv(outfilepathToSave)
            
        except ValueError:
            print "No reads in %s"%bedFile

print "Generate sorted gene lists by gene type."
logOpen.write("Generate sorted gene lists by gene type.\n")

bedFile_pc=outfilepath+"clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed"
geneCounts_pc=sortFilteredBed(bedFile_pc) 
outfilepathToSave=outfilepath + '/PlotData_ReadsPerGene_proteinCoding'
geneCounts_pc.to_csv(outfilepathToSave)

bedFile_linc=outfilepath+"clipGenes_lincRNA_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed"
geneCounts_linc=sortFilteredBed(bedFile_linc)
outfilepathToSave=outfilepath + '/PlotData_ReadsPerGene_lincRNA'
geneCounts_linc.to_csv(outfilepathToSave)

CLIPPERlowFDRcenters=getBedCenterPoints(CLIPPERlowFDR)
allLowFDRCentersBedGraph=makeBedGraph(CLIPPERlowFDRcenters,genomeFile)
bedFile_sno=getSnoRNAreads(CLIPPERlowFDRcenters,snoRNAindex)
geneCounts_sno=countSnoRNAs(bedFile_sno) 
outfilepathToSave=outfilepath + '/PlotData_ReadsPerGene_snoRNA'
geneCounts_sno.to_csv(outfilepathToSave)
    
remaining=[f for f in glob.glob(outfilepath+"*_LowFDRreads.bed") if 'lincRNA' not in f and 'proteinCoding' not in f and 'snoRNA' not in f]
countRemainingGeneTypes(remaining)

# <codecell>

def makeClusterCenter(windowsFile):
    # Usage: Generate a file of cluster centers.
    # Input: Raw CLIPper output file.
    # Output: File with coordinates for the center of each CLIPper cluster.
    cleanBed = cleanBedFile(windowsFile)
    centers=cleanBed.replace('.bed','.clusterCenter')
    f = open(centers, 'w')
    with open(cleanBed, 'r') as infile:
        for line in infile:
            elementList = line.strip().split('\t')
            diff=abs(int((int(elementList[1])-int(elementList[2]))/2))
            f.write(elementList[0]+'\t'+str(int(elementList[1])+diff)+'\t'+str(int(elementList[1])+diff+1)+'\n')
    f.close()
    return centers

def getClusterIntensity(bedGraph,centerCoordinates):
    # Usage: Generate a matrix of read itensity values around CLIPper cluster center.
    # Input: BedGraph and cluster center file.
    # Output: Generates a matrix, which is passed into R.
    program=os.getcwd() + '/bin/grep_chip-seq_intensity.pl'
    program2='wait'
    proc=subprocess.Popen(['perl',program, centerCoordinates, bedGraph],)
    proc.communicate()
    logOpen.write("Waiting for Cluster Intensity file completion...\n")
    proc2=subprocess.Popen(program2,shell=True)
    proc2.communicate()
    
print "Get binding intensity around cluster centers."
logOpen.write("Get binding intensity around cluster centers.\n")
bedGraphCLIPin=makeBedGraph(CLIPPERin,genomeFile)
centerCoordinates=makeClusterCenter(CLIPperOutBed) 
getClusterIntensity(bedGraphCLIPin,centerCoordinates)

# <codecell>

def partitionReadsByUTR(infile,UTRmask,utrReads,notutrReads):
    program = 'intersectBed'
    outfh = open(utrReads,'w')
    proc = subprocess.Popen([program,'-a',infile,'-b',UTRmask,'-u','-s'],stdout=outfh)
    proc.communicate()
    outfh.close()
    outfh = open(notutrReads,'w')
    proc = subprocess.Popen([program,'-a',infile,'-b',UTRmask,'-v','-s'],stdout=outfh)
    proc.communicate()
    outfh.close()

def extractUTRs(bedIn,fivePUTRBed,threePUTRBed,cdsBed):
    # Usage: Extract all UTR specific reads from the input file.
    # Input: .bed file
    # Output: Mutually exclusive partitions of the input file.
    fivePreads = bedIn.replace('.bed', '_5p.bed')
    notFivePreads = bedIn.replace('.bed', '_NOT5p.bed')
    partitionReadsByUTR(bedIn,fivePUTRBed,fivePreads,notFivePreads)
    threePreads = bedIn.replace('.bed', '_3p.bed')
    notThreePreads = bedIn.replace('.bed', '_NOT3p.bed')
    partitionReadsByUTR(notFivePreads,threePUTRBed,threePreads,notThreePreads)
    CDSreads = bedIn.replace('.bed', '_cds.bed')
    notCDSreads = bedIn.replace('.bed', '_NOTcds.bed')
    partitionReadsByUTR(notThreePreads,cdsBed,CDSreads,notCDSreads)
    return (fivePreads,notFivePreads,CDSreads,notCDSreads,threePreads,notThreePreads)

print "Intron and UTR analysis."
logOpen.write("Intron and UTR analysis.\n")
fivePreads,notFivePreads,CDSreads,notCDSreads,threePreads,notThreePreads=extractUTRs(filteredProteinCodingCenters,fivePUTRBed,threePUTRBed,cdsBed)
geneCounts_5p=sortFilteredBed(fivePreads) 
geneCounts_3p=sortFilteredBed(threePreads) 
geneCounts_cds=sortFilteredBed(CDSreads) 

outfilepathToSave=outfilepath+'/PlotData_ReadsPerGene_5pUTR'
geneCounts_5p.to_csv(outfilepathToSave)
outfilepathToSave=outfilepath+'/PlotData_ReadsPerGene_3pUTR'
geneCounts_3p.to_csv(outfilepathToSave)
outfilepathToSave=outfilepath+'/PlotData_ReadsPerGene_CDS'
geneCounts_cds.to_csv(outfilepathToSave) 

# <codecell>

def makeTab(bedGraph,genesFile,sizesFile):
    program = os.getcwd() + '/bin/bedGraph2tab.pl'
    program2 = 'wait'
    outfile=bedGraph.replace('.bedgraph','.tab')
    proc = subprocess.Popen(['perl',program,genesFile,sizesFile,bedGraph,outfile],)
    proc.communicate()
    proc2 = subprocess.Popen(program2,shell=True)
    proc2.communicate()
    return outfile

def makeAvgGraph(bedGraph,utrFile,genesFile,sizesFile):
    # Usage: Generate a matrix of read itensity values across gene body.
    # Input: BedGraph.
    # Output: Generates two matricies.
    program= os.getcwd() + '/bin/averageGraph_scaled_tab.pl'
    program2 = 'wait'
    tabFile=makeTab(bedGraph,genesFile,sizesFile)
    outhandle=tabFile.replace('.tab','_UTRs')
    proc = subprocess.Popen(['perl',program,utrFile,tabFile,tabFile,outhandle],)
    proc.communicate()
    proc2 = subprocess.Popen(program2,shell=True)
    proc2.communicate()

print "Gene body analysis."
logOpen.write("Gene body analysis.\n")
bedGraphProtein=makeBedGraph(bedFile_pc,genomeFile)
makeAvgGraph(bedGraphProtein,utrFile,genesFile,sizesFile)

# <codecell>

def getGeneStartStop(bedFile,geneRef):
    try:
        bf=pd.DataFrame(pd.read_table(bedFile,header=None))
        bf.columns=['Chr','Start','End','ReadName','Q','Strand','CLIPper_winChr','CLIPper_winStart','CLIPper_winEmd','CLIPper_winaName','CLIPper_winP','CLIPper_winStrand']
        bf['geneName']=bf['CLIPper_winaName'].apply(lambda x: x.split('_')[0])
        merge=pd.merge(geneRef,bf,left_on='Ensembl Gene ID',right_on='geneName')
        ncRNA_startStop=merge[['Ensembl Gene ID','Gene Start (bp)','Gene End (bp)','Start','End','Strand']]
        outfilepathToSave=bedFile.replace(".bed",".geneStartStop")
        ncRNA_startStop.to_csv(outfilepathToSave)
    except ValueError:
        print "No reads in %s"%bedFile

print "ncRNA gene body anaysis."
geneStartStopRepo=os.getcwd()+'/docs/all_genes.txt'
geneRef=pd.DataFrame(pd.read_table(geneStartStopRepo))
remaining=[f for f in glob.glob(outfilepath+"*_LowFDRreads.bed") if 'lincRNA' not in f and 'proteinCoding' not in f and 'snoRNA' not in f]
for bedFile in remaining:
    st_stop=getGeneStartStop(bedFile,geneRef)

# lincRNA file processing
bedFile_linc=outfilepath+"clipGenes_lincRNA_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed"
bf=pd.DataFrame(pd.read_table(bedFile_linc,header=None))
bf.columns=['Chr','Start','Stop','CLIPper_name','Q','Strand']
bf['geneName']=bf['CLIPper_name'].apply(lambda x: x.split('_')[0])
merge=pd.merge(geneRef,bf,left_on='Ensembl Gene ID',right_on='geneName')
ncRNA_startStop=merge[['Ensembl Gene ID','Gene Start (bp)','Gene End (bp)','Start','Stop','Strand']]
outfilepathToSave=bedFile_linc.replace(".bed",".geneStartStop")
ncRNA_startStop.to_csv(outfilepathToSave)

# <codecell>

def makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation):
    repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
    repeat_genome_bases=repeat_genome[1]
    repeat_genome_size=len(repeat_genome[1])
    repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
    repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
    repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 # Python list extraction is not end index inclusive; to extract sequence, use end + 1. 
    return (repeat_genome_bases,repeatAnnotDF)

def readBed(path):
    bedFile = pd.read_table(path,dtype=str,header=None)
    bedFile.columns=['Index','Start','Stop','Name','QS','Strand']
    bedFile['Start']=bedFile['Start'].astype(int)
    return bedFile

print "Record repeat RNA."
repeat_genome_bases,repeatAnnotDF=makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation)
repeatAnnotDF.set_index('Name',inplace=True,drop=False)
# Get merged data for repeat index.
repeatMerged=glob.glob(outfilepath+"*repeat_allreads.mergedRT.bed")
rep=pd.read_table(repeatMerged[0],dtype=str,header=None)
rep.columns=['Rep_index','Start','Stop','Read_name','Q','Strand']
rep['RT_stop']=rep['Start'].astype(int)+expand
for ix in repeatAnnotDF.index:
    end=repeatAnnotDF.loc[ix,'IndexEnd']
    repName=repeatAnnotDF.loc[ix,'Name']
    gene_hits=rep[(rep['RT_stop']<int(repeatAnnotDF.loc[ix,'IndexEnd']))&(rep['RT_stop']>int(repeatAnnotDF.loc[ix,'IndexStart']))]
    gene_hits['Repeat_End']=repeatAnnotDF.loc[ix,'IndexEnd']
    gene_hits['Repeat_Start']=repeatAnnotDF.loc[ix,'IndexStart']
    outfilepathToSave=outfilepath + '/PlotData_RepeatRNAreads_%s'%repName
    gene_hits.to_csv(outfilepathToSave)

# <codecell>

def makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation):
    repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
    repeat_genome_bases=repeat_genome[1]
    repeat_genome_size=len(repeat_genome[1])
    repeatAnnotDF=pd.DataFrame(pd.read_table(repeatAnnotation,header=None))
    repeatAnnotDF.columns=['Name','Length','IndexStart','IndexEnd']
    repeatAnnotDF['End_for_extraction']=repeatAnnotDF['IndexEnd']+1 # Python list extraction is not end index inclusive; to extract sequence, use end + 1. 
    return (repeat_genome_bases,repeatAnnotDF)

repeat_genome_bases,repeatAnnotDF=makeRepeatAnnotation(repeatGenomeBuild,repeatAnnotation)

# <codecell>

def lineCount(filename):
	i=0
	with open(filename) as f:
		for i,l in enumerate(f):
			pass
	return i+1

def plot_ReadAccounting(outfilepath,sampleName):
    rawRead1=infilepath+sampleName+'_R1.fastq'
    rawRead2=infilepath+sampleName+'_R2.fastq'
    reads3pTrim=[outfilepath+sampleName+'_R1_3ptrimmed.fastq',outfilepath+sampleName+'_R2_3ptrimmed.fastq']
    readsFilter=[outfilepath+sampleName+'_R1_3ptrimmed_filter.fastq',outfilepath+sampleName+'_R2_3ptrimmed_filter.fastq']
    readsNoDupes=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe.fastq',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe.fastq']
    readsMappedReapeat=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_mappedTorepeat_withDupes.bed',outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_mappedTorepeat_withDupes.bed']
    readsMappedHg19=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes.bed'%index_tag,outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes.bed'%index_tag]
    readsMappedBlacklist=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes.bed'%index_tag,outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes.bed'%index_tag]
    readsMappedRepeatMask=[outfilepath+sampleName+'_R1_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes_noBlacklist_noRepeat.bed'%index_tag,outfilepath+sampleName+'_R2_3ptrimmed_filter_nodupe_5ptrimmed_notMappedTorepeat_mappedTo%s_withDupes_noBlacklist_noRepeat.bed'%index_tag]
    clipperIN=outfilepath+sampleName+'_threshold=%s_%s_allreads.mergedRT_CLIPPERin.bed'%(threshold,index_tag)
    clipperOUT=outfilepath+sampleName+'_threshold=%s_%s_allreads.mergedRT_CLIP_clusters_lowFDRreads.bed'%(threshold,index_tag)
    fileNames=['Raw (R1)','Raw (R2)','3p Trim (R1)','3p Trim (R2)','Filter (R1)','Filter (R2)','No dupes (R1)','No dupes (R2)','RepeatMapped (R1)','RepeatMapped (R2)','Hg19Mapped (R1)','Hg19Mapped (R2)','Blacklist (R1)','Blacklist (R2)','RepeatMask (R1)','RepeatMask (R2)','ClipperIn','ClipperOut']
    filesToCount=[rawRead1,rawRead2,reads3pTrim[0],reads3pTrim[1],readsFilter[0],readsFilter[1],readsNoDupes[0],readsNoDupes[1],readsMappedReapeat[0],readsMappedReapeat[1],readsMappedHg19[0],readsMappedHg19[1],readsMappedBlacklist[0],readsMappedBlacklist[1],readsMappedRepeatMask[0],readsMappedRepeatMask[1],clipperIN,clipperOUT]
    
    counts=[]
    counter=0
    for fileString in filesToCount:
        temp=lineCount(fileString)
        if counter < 8:
            temp=temp/4 # Fastq files
        counts=counts+[temp]
        counter += 1

    ind = np.arange(len(counts)) + 0.5
    plt.barh(ind,list(reversed(np.log10(np.array(counts)))),align='center',color='blue')
    plt.xlabel('log10(Counts per file)',fontsize=5)
    locs,pltlabels = plt.xticks(fontsize=5)
    plt.setp(pltlabels, rotation=90, fontsize=5)
    plt.yticks(ind,list(reversed(fileNames)),fontsize=5)
    plt.tick_params(axis='yticks',labelsize=5) 
    ax=plt.gca()
    for line in ax.get_yticklines():
        line.set_markersize(0)
    plt.title('Read counts',fontsize=5)
    
    readDF=pd.DataFrame()
    readDF['File_name']=fileNames
    readDF['Reads_per_file']=counts
    outfilepathToSave=outfilepath + '/PlotData_ReadsPerPipeFile'
    readDF.to_csv(outfilepathToSave)

plt.subplot(2,3,1) 
plot_ReadAccounting(outfilepath,sampleName)

# <codecell>

def plot_BoundGeneTypes(outfilepath,sampleName):
    record=pd.DataFrame()   
    # Exclude specific files (e.g., UTR-specific reads).
    geneListToPlot=[f for f in glob.glob(outfilepath+'PlotData_ReadsPerGene_*') if '5pUTR' not in f and '3pUTR' not in f and 'CDS' not in f]
    for boundGenes in geneListToPlot:
        glist=pd.read_csv(boundGenes,header=None)
        glist.columns=['GeneName','Count']
        gName=boundGenes.split('_')[-1]
        record.loc[gName,'genesBound']=glist.shape[0]
        record.loc[gName,'totalReads']=glist['Count'].sum()
    record.sort('genesBound',inplace=True)
    outfilepathToSave=outfilepath + '/PlotData_ReadAndGeneCountsPerGenetype'
    record.to_csv(outfilepathToSave)
    ind = np.arange(record.shape[0]) + 0.5
    plt.bar(ind,record['genesBound'],align='center',color='blue')
    locs,pltlabels = plt.yticks(fontsize=5)
    locs,pltlabels = plt.xticks(ind,record.index,fontsize=5)
    plt.setp(pltlabels, rotation=90, fontsize=5)
    plt.tick_params(axis='xticks',labelsize=5) 
    ax=plt.gca()
    for line in ax.get_xticklines():
        line.set_markersize(0)
    plt.ylabel('Number of genes bound',fontsize=5)
    plt.tick_params(axis='yticks',labelsize=5)
    plt.title('Bound genes by class',fontsize=5)
    
plt.subplot(2,3,6)
plot_BoundGeneTypes(outfilepath,sampleName)

# <codecell>

def plot_ReadsPerCluster(outfilepath,sampleName):
    readPerCluster=outfilepath+sampleName+'_threshold=%s_%s_allreads.mergedRT_CLIP_clusters.readsPerCluster'%(threshold,index_tag)
    clust=pd.DataFrame(pd.read_table(readPerCluster,header=None))
    clust.columns=['ReadsPerCluster']
    clust=clust['ReadsPerCluster']
    interval=10
    bins=range(min(clust)-10,max(clust)+10,interval)
    hist,bins=np.histogram(clust,bins=bins)
    width=0.7*(bins[1]-bins[0]) 
    center=(bins[:-1] + bins[1:])/2 
    plt.bar(center, hist,align='center',width=width)
    locs,pltlabels = plt.yticks(fontsize=5)
    locs,pltlabels = plt.xticks(center,center,fontsize=5)
    plt.setp(pltlabels, rotation=90, fontsize=3.5)
    plt.tick_params(axis='yticks',labelsize=5) 
    plt.xlabel('Reads per cluster (bin=%s)'%interval,fontsize=5)
    plt.ylabel('Frequency (RT stop count)',fontsize=5)
    plt.title('Reads per cluster',fontsize=5)
    plt.xlim(0,100) # Make the histogram easy to view.
    # plt.xlim(-interval,np.max(center)+interval)
    
plt.subplot(2,3,2)
plot_ReadsPerCluster(outfilepath,sampleName)

# <codecell>

def plot_ClusterSizes(outfilepath,sampleName):
    clipClusters=outfilepath+sampleName+'_threshold=%s_%s_allreads.mergedRT_CLIP_clusters'%(threshold,index_tag)
    clust=pd.DataFrame(pd.read_table(clipClusters,header=None,skiprows=1))
    clust.columns=['chr','start','end','name','score','strand','m1','m2']
    clust['clusterSize']=clust['start']-clust['end']
    clust['clusterSize']=clust['clusterSize'].apply(lambda x: math.fabs(x))
    plt.boxplot(clust['clusterSize'])
    plt.tick_params(axis='x',labelbottom='off') 
    ax=plt.gca()
    for line in ax.get_xticklines():
        line.set_markersize(0)
    plt.ylabel('Cluster length (bases)',fontsize=5)
    locs,pltlabels = plt.yticks(fontsize=5)
    plt.title('Cluster size',fontsize=5)

plt.subplot(2,3,3)
plot_ClusterSizes(outfilepath,sampleName)

# <codecell>

def plot_clusterBindingIntensity(outfilepath,sampleName):
    clusterCenterHeatmap=outfilepath+sampleName+'_threshold=%s_%s_allreads.mergedRT_CLIP_clusters_cleaned_sorted.clusterCenter_heatmap.txt'%(threshold,index_tag)
    hmap=pd.DataFrame(pd.read_table(clusterCenterHeatmap,header=None,skiprows=1))
    hmap_vals=hmap.ix[:,1:]
    sums=hmap_vals.sum(axis=1)
    hmap_vals=hmap_vals.loc[np.argsort(sums),:]
    plt.ylim(0,hmap_vals.shape[0])
    p=plt.pcolormesh(np.array(hmap_vals),cmap='Blues')
    plt.tick_params(axis='x',labelbottom='off') 
    plt.xlabel('Cluster position',fontsize=5)
    locs,pltlabels = plt.yticks(fontsize=5)
    plt.ylabel('Cluster number',fontsize=5)
    plt.title('Read distribution',fontsize=5)

plt.subplot(2,3,4)
plot_clusterBindingIntensity(outfilepath,sampleName)

# <codecell>

def readUTRfile(path):
    geneCounts=pd.read_csv(path,header=None)
    geneCounts.columns=['Gene_name','Count']
    return geneCounts

def plot_readsBymRNAregion(outfilepath,sampleName): 
    pc_5pReads=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_5pUTR')['Count'].sum()
    pc_3pReads=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_3pUTR')['Count'].sum()
    pc_CDSReads=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_CDS')['Count'].sum()
    non_intronic=pc_5pReads+pc_3pReads+pc_CDSReads
    allProteinCoding=outfilepath +'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved.bed'
    all_pc=pd.DataFrame(pd.read_table(allProteinCoding,header=None))
    pc_allReads=all_pc.shape[0]
    v=[float(pc_allReads-non_intronic)/pc_allReads,float(pc_5pReads)/pc_allReads,float(pc_CDSReads)/pc_allReads,float(pc_3pReads)/pc_allReads]
    pie_wedges=ax.pie(v,labels=["Intronic","5p UTR","CDS","3pUTR"],labeldistance=1.1,autopct='%1.1f%%')
    plt.rcParams['font.size']=5
    for wedge in pie_wedges[0]:
        wedge.set_edgecolor('black')
        wedge.set_lw(1)

ax=plt.subplot(2,3,5)
plot_readsBymRNAregion(outfilepath,sampleName)

# <codecell>

fig1=plt.figure(1)

plt.subplot(2,3,1) 
plot_ReadAccounting(outfilepath,sampleName)
plt.subplot(2,3,2)
plot_ReadsPerCluster(outfilepath,sampleName)
plt.subplot(2,3,3)
plot_ClusterSizes(outfilepath,sampleName)
plt.subplot(2,3,4)
plot_clusterBindingIntensity(outfilepath,sampleName)
ax=plt.subplot(2,3,5)
plot_readsBymRNAregion(outfilepath,sampleName)
plt.subplot(2,3,6)
plot_BoundGeneTypes(outfilepath,sampleName)
fig1.tight_layout()

fig1.savefig(outfilepath+'Figure1.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig1.savefig(outfilepath+'Figure1.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

def plot_mRNAgeneBodyDist(outfilepath,sampleName):
    averageGraph=outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved_cleaned_sorted_UTRs_scaled_cds200_abt0_averageGraph.txt'
    hmap=pd.DataFrame(pd.read_table(averageGraph,header=None,skiprows=1))
    hmap=hmap.set_index(0)
    avgTrace=hmap.loc['treat',:]
    plt.plot(avgTrace,color='blue',linewidth='2')
    plt.vlines(200,0,np.max(avgTrace),linestyles='dashed')
    plt.vlines(400,0,np.max(avgTrace),linestyles='dashed')
    plt.ylim(0,np.max(avgTrace))
    plt.tick_params(axis='x',labelbottom='off') 
    plt.xlabel('mRNA gene body (5pUTR, CDS, 3pUTR)')
    plt.ylabel('Read density')
    plt.tick_params(axis='y',labelsize=5) 
    plt.title('CLIP signal across average mRNA transcript.',fontsize=5)

plt.subplot2grid((2,3),(0,0),colspan=3)
plot_mRNAgeneBodyDist(outfilepath,sampleName)

# <codecell>

def convertENBLids(enst_name):
    ensg_name=ensemblGeneAnnot.loc[enst_name,'name2']
    return ensg_name

def getUTRbindingProfile(utr,hmap_m):
    if utr=='5p':
        ix=(hmap_m[range(201,601)].sum(axis=1)==0)&(hmap_m[range(1,201)].sum(axis=1)>0)
        screen=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_5pUTR')
    elif utr=='3p':
        ix=(hmap_m[range(1,401)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)>0)
        screen=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_3pUTR')
    else:
        ix=(hmap_m[range(1,201)].sum(axis=1)==0)&(hmap_m[range(401,601)].sum(axis=1)==0)&(hmap_m[range(201,401)].sum(axis=1)>0)
        screen=readUTRfile(outfilepath+'/PlotData_ReadsPerGene_CDS')
        
    # Ensure all genes are also identified in pre-allocated gene lists.
    hmap_m_utrSpec=hmap_m.ix[ix,:]
    hmap_m_utrSpec_filter=pd.merge(hmap_m_utrSpec,screen,left_on='ENSG_ID',right_on='Gene_name',how='inner')
    sums=hmap_m_utrSpec_filter[range(1,601)].sum(axis=1)
    hmap_m_utrSpec_filter=hmap_m_utrSpec_filter.loc[np.argsort(sums),:]
    return hmap_m_utrSpec_filter

def plot_geneBodyPartition(outfilepath,sampleName):
    treatMatrix=outfilepath+'clipGenes_proteinCoding_LowFDRreads_centerCoord_snoRNAremoved_miRNAremoved_cleaned_sorted_UTRs_scaled_cds200_abt0_treatmatrix.txt'
    hmap=pd.DataFrame(pd.read_table(treatMatrix,header=None,skiprows=1))
    
    # Ensure genes recoverd from this analysis are indepdently indentified using partitioning of CLIPper cluster data.
    hmap['ENSG_ID']=hmap.ix[:,0].apply(convertENBLids)
    bound_pc = outfilepath+'clipGenes_proteinCoding'
    pc_genes=pd.DataFrame(pd.read_table(bound_pc,header=None,))
    pc_genes.columns=['ENSG_ID']
    hmap_m=pd.merge(hmap,pc_genes,left_on='ENSG_ID',right_on='ENSG_ID',how='inner') 
    
    # Isolate intronic bound genes.
    tosave=outfilepath+'PlotData_ExclusiveBound_Intronic' 
    intronicBoundGenes=list(set(pc_genes['ENSG_ID'])-set(hmap_m['ENSG_ID']))
    np.savetxt(tosave,np.array(intronicBoundGenes),fmt="%s")
    
    # UTR specific genes.
    geneTypes=['5p','cds','3p'] 
    depth=50
    for i in range(0,3):    
        utrMatrix=getUTRbindingProfile(geneTypes[i],hmap_m)
        tosave=outfilepath+'PlotData_ExclusiveBound_%s'%geneTypes[i] 
        np.savetxt(tosave,utrMatrix['ENSG_ID'],fmt="%s")
        plt.subplot2grid((2,3),(1,i),colspan=1)
        dataToPlot=utrMatrix[range(1,601)]
        p=plt.pcolormesh(np.array(dataToPlot)[-depth:-1,:],cmap='Blues')
        plt.title(geneTypes[i],fontsize=5)
        plt.vlines(200,0,depth,linestyles='dashed')
        plt.vlines(400,0,depth,linestyles='dashed')
        plt.tick_params(axis='x',labelbottom='off') 
        plt.tick_params(axis='y',labelleft='off') 
        plt.ylim(0,depth)
        plt.ylabel('Ranked genes (highest on bottom)',fontsize=5)
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        plt.title('%s specific genes: %s'%(geneTypes[i],np.unique(utrMatrix['ENSG_ID']).shape[0]),fontsize=5)
 
ensemblGeneAnnot=pd.DataFrame(pd.read_table(genesFile))
ensemblGeneAnnot=ensemblGeneAnnot.set_index('name') # Make ENST the index
plot_geneBodyPartition(outfilepath,sampleName)

# <codecell>

fig2=plt.figure(2)
plt.subplot2grid((2,3),(0,0),colspan=3)
plot_mRNAgeneBodyDist(outfilepath,sampleName)
plot_geneBodyPartition(outfilepath,sampleName)
fig2.tight_layout()
fig2.savefig(outfilepath+'Figure2.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig2.savefig(outfilepath+'Figure2.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

def plot_repeatRNA(outfilepath,sampleName):
    repeat_genome=np.genfromtxt(repeatGenomeBuild,dtype='string')
    repeat_genome_bases=repeat_genome[1]
    
    repFiles=glob.glob(outfilepath + '/PlotData_RepeatRNAreads_*')
    repFiles=[repFile for repFile in repFiles if 'rDNA' not in repFile]
    
    plotDim=math.ceil(math.sqrt(len(repFiles)))
    i=0
    for path in repFiles:
        name=path.split('RepeatRNAreads_')[-1]
        try:
            # Read in each RT stop file
            hits_per_rep=pd.read_csv(path)
            RTpositions=hits_per_rep['RT_stop']
            start=hits_per_rep.loc[0,'Repeat_Start']
            end=hits_per_rep.loc[0,'Repeat_End']
            # Histogram of RT stops across gene body
            bins=range(start,end+2,1)
            hist,bins=np.histogram(RTpositions,bins=bins)
            width=0.7*(bins[1]-bins[0])
            center=(bins[:-1] + bins[1:])/2
            # Normalize
            histPlot=np.array(hist,dtype=float)
            histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
            # Subplot
            plt.subplot(plotDim,plotDim,i+1)
            plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
            plt.tick_params(axis='x',labelsize=2.5) 
            plt.tick_params(axis='y',labelsize=2.5)  
            plt.title('RT stops for %s: %s'%(name,len(RTpositions)),fontsize=5)
            plt.xlim(start,end)  
            # Record data
            storageDF=pd.DataFrame()
            sequence=repeat_genome_bases[start:end+1]
            storageDF['Sequence']=pd.Series(list(sequence))
            readsPerBase=np.array(list(hist))
            readsPerBaseNorm=np.array(list(histPlot))
            storageDF['RT_stops']=readsPerBase
            storageDF['RT_stops_norm']=readsPerBaseNorm              
            outfilepathToSave=outfilepath +'/PlotData_RepeatRNAHist_%s'%name
            storageDF.to_csv(outfilepathToSave)
            i+=1
        except:
            print "No reads for repeatRNA %s"%name              
    plt.tight_layout()

fig3=plt.figure(3)
plot_repeatRNA(outfilepath,sampleName)
fig3.tight_layout()
fig3.savefig(outfilepath+'Figure3.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig3.savefig(outfilepath+'Figure3.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

def plot_rDNA(outfilepath,sampleName):
    plt.subplot2grid((3,3),(0,0),colspan=3)
    name='rDNA'
    rDNA=glob.glob(outfilepath + 'PlotData_RepeatRNAreads_rDNA')
    hits_per_rep=pd.read_csv(rDNA[0])
    RTpositions=hits_per_rep['RT_stop']
    start=hits_per_rep.loc[0,'Repeat_Start']
    end=hits_per_rep.loc[0,'Repeat_End']
    bins=range(start,end+2,1)
    hist,bins=np.histogram(RTpositions,bins=bins)
    width=0.7*(bins[1]-bins[0])
    center=(bins[:-1]+bins[1:])/2
    histPlot=np.array(hist,dtype=float)
    histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
    plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
    plt.tick_params(axis='x',labelsize=2.5) 
    plt.tick_params(axis='y',labelsize=2.5)  
    plt.title('RT stops for %s: %s'%(name,len(RTpositions)),fontsize=5)
    plt.xlim(start,end)  
    # Record data
    storageDF=pd.DataFrame()
    sequence=repeat_genome_bases[start:end+1]
    storageDF['Sequence']=pd.Series(list(sequence))
    readsPerBase=np.array(list(hist))
    readsPerBaseNorm=np.array(list(histPlot))
    storageDF['RT_stops']=readsPerBase
    storageDF['RT_stops_norm']=readsPerBaseNorm              
    outfilepathToSave=outfilepath +'/PlotData_RepeatRNAHist_%s'%name
    storageDF.to_csv(outfilepathToSave)

    # Features of rDNA with respect to start of the bowtie index (index=0)
    rRNAstart=start
    plt.axvspan(start18s+rRNAstart,end18s+rRNAstart,facecolor='g',alpha=0.5)
    plt.axvspan(start5s+rRNAstart,end5s+rRNAstart,facecolor='r',alpha=0.5)
    plt.axvspan(start28s+rRNAstart,end28s+rRNAstart,facecolor='b',alpha=0.5)
    
    # Generate histogram for transcribed region
    plt.subplot2grid((3,3),(1,0),colspan=3)
    datarDNAOnly=RTpositions-start
    bins=range((start-start),(end-start+2),1)
    hist,bins=np.histogram(datarDNAOnly,bins=bins)
    width=0.7*(bins[1]-bins[0])
    center=(bins[:-1] + bins[1:])/2
    histPlot=np.array(hist,dtype=float)
    histPlot=np.array(histPlot/float(len(RTpositions)),dtype=float)
    plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.45)
    plt.tick_params(axis='x',labelsize=2.5) 
    plt.tick_params(axis='y',labelsize=2.5)  
    plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
    plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
    plt.axvspan(start18s,end18s,facecolor='g',alpha=0.5)
    plt.axvspan(start5s,end5s,facecolor='r',alpha=0.5)
    plt.axvspan(start28s,end28s,facecolor='b',alpha=0.5)
    plt.xlim(0,rRNAend)
    
    # Individual regions 
    plt.subplot2grid((3,3),(2,0),colspan=1)
    plt.bar(center,histPlot,align='center',width=width,color='green',alpha=0.75)
    plt.xlim(start18s,end18s)
    plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
    plt.ylabel('Normalized RT stop / bin',fontsize=2.5)
    plt.tick_params(axis='x',labelsize=5) 
    plt.tick_params(axis='y',labelsize=5) 
    plt.title('18s Region',fontsize=5)
    plt.subplot2grid((3,3),(2,1),colspan=1)
    plt.bar(center,histPlot,align='center',width=width,color='red',alpha=0.75)
    plt.xlim(start5s,end5s)
    plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
    plt.tick_params(axis='x',labelsize=5) 
    plt.tick_params(axis='y',labelsize=5) 
    plt.title('5.8s Region',fontsize=5)
    plt.subplot2grid((3,3),(2,2),colspan=1)
    plt.bar(center,histPlot,align='center',width=width,color='blue',alpha=0.75)
    plt.xlim(start28s,end28s)
    plt.xlabel('rRNA locus position (bin=1 base)',fontsize=5)
    plt.tick_params(axis='x',labelsize=5) 
    plt.tick_params(axis='y',labelsize=5)  
    plt.title('28s Region',fontsize=5)
    plt.tight_layout()

fig4=plt.figure(4)
plot_rDNA(outfilepath,sampleName)
fig4.tight_layout()
fig4.savefig(outfilepath+'Figure4.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig4.savefig(outfilepath+'Figure4.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

def getBindingFrac(type_specific):
    # 5' position on the negative strand is snoRNA stop coordinate.
    neg_data=type_specific[type_specific['strand_snoRNA']=='-']
    neg_data['diff']=np.abs(neg_data['Stop_snoRNA']-neg_data['Start']) 
    neg_data['frac']=neg_data['diff']/(neg_data['Stop_snoRNA']-neg_data['Start_snoRNA'])
    # 5' position on the positive strand is snoRNA start coordinate.
    pos_data=type_specific[type_specific['strand_snoRNA']=='+']
    pos_data['diff']=np.abs(pos_data['Start_snoRNA']-pos_data['Start'])
    pos_data['frac']=pos_data['diff']/(pos_data['Stop_snoRNA']-pos_data['Start_snoRNA'])
    DF_snoProfile=pd.concat([neg_data,pos_data])
    return DF_snoProfile

print "snoRNA gene body anaysis."
# logOpen.write("Gene body analysis.\n")
bf_sno=pd.read_table(outfilepath+"clipGenes_snoRNA_LowFDRreads.bed",header=None)
bf_sno.columns=['Chr','Start','End','CLIPper_name','Q','Strand','Chr_snoRNA','Start_snoRNA','Stop_snoRNA','name_snoRNA','Type','strand_snoRNA']
snoTypes=pd.DataFrame(bf_sno.groupby('Type').size())
snoTypes.columns=['Reads']
snoTypes['Fraction']=snoTypes['Reads']/snoTypes['Reads'].sum(axis=1)
outfilepathToSave=outfilepath+'/PlotData_readsPerSnoRNAType'
snoTypes.to_csv(outfilepathToSave)

snoTypesAndGenes=pd.DataFrame(bf_sno.groupby(['Type','name_snoRNA']).size())
snoTypesAndGenes.columns=['Count_per_gene']
outfilepathToSave=outfilepath+'/PlotData_geneStatsPerSnoRNAType'
snoTypesAndGenes.to_csv(outfilepathToSave)

fig5=plt.figure(5)
ax=plt.subplot(2,2,1)
pie_wedges=ax.pie(snoTypes['Fraction'],labels=snoTypes.index,labeldistance=1.1,autopct='%1.1f%%')
plt.rcParams['font.size']=5
for wedge in pie_wedges[0]:
    wedge.set_edgecolor('black')
    wedge.set_lw(1)

i=2
for sType in set(bf_sno['Type']):
    type_specific=bf_sno[bf_sno['Type']==sType]
    sno_profile=getBindingFrac(type_specific)
    
    if sType=='C':
        title="C/D_box"
    elif sType=='H':
        title="H/ACA_box"
    else:
        title="scaRNA"
    
    outfilepathToSave=outfilepath+'/PlotData_snoRNAReadDist_%s'%sType
    sno_profile.to_csv(outfilepathToSave)
    
    plt.subplot(2,2,i)
    bins=np.arange(0,1,0.01)
    hist,bins=np.histogram(sno_profile['frac'],bins=bins)
    hist=np.array(hist/float(sno_profile['frac'].shape[0]),dtype=float)
    width=0.7*(bins[1]-bins[0])
    center=(bins[:-1] + bins[1:])/2
    plt.bar(center,hist,align='center',width=width,color='blue',alpha=0.75)
    plt.tick_params(axis='x',labelsize=5) 
    plt.tick_params(axis='y',labelsize=5)  
    plt.xlabel('Fraction of gene body (5p - 3p)',fontsize=5)
    plt.title('Binding profile for %s'%title,fontsize=5)
    plt.xlim([0,1])
    
    # Record data
    storageDF=pd.DataFrame()
    storageDF['bins']=pd.Series(bins)
    storageDF['hist']=pd.Series(hist)
    outfilepathToSave=outfilepath+'/PlotData_snoRNAhistogram_%s'%sType
    storageDF.to_csv(outfilepathToSave)
    
    i+=1

fig5.tight_layout()
fig5.savefig(outfilepath+'Figure5.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig5.savefig(outfilepath+'Figure5.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

def getncRNABindingFrac(type_specific):
    # 5' position on the negative strand is snoRNA stop coordinate.
    neg_data=type_specific[type_specific['Strand']=='-']
    neg_data['diff']=np.abs(neg_data['Gene End (bp)']-neg_data['RT_stop']) 
    neg_data['frac']=neg_data['diff']/(neg_data['Gene End (bp)']-neg_data['Gene Start (bp)'])
    # 5' position on the positive strand is snoRNA start coordinate.
    pos_data=type_specific[type_specific['Strand']=='+']
    pos_data['diff']=np.abs(pos_data['Gene Start (bp)']-pos_data['RT_stop'])
    pos_data['frac']=pos_data['diff']/(pos_data['Gene End (bp)']-pos_data['Gene Start (bp)'])
    DF_ncRNAProfile=pd.concat([neg_data,pos_data])
    return DF_ncRNAProfile

print "ncRNA gene body anaysis."
st_stopFiles=glob.glob(outfilepath+"*.geneStartStop")
st_stopFiles=[f for f in st_stopFiles if 'rRNA' not in f]
fig6=plt.figure(6)
plotDim=math.ceil(math.sqrt(len(st_stopFiles)))
i=1
for st_file in st_stopFiles:
    name=st_file.split('clipGenes_')[1].split('_LowFDRreads')[0]
    tmp=pd.read_csv(st_file)
    tmp['RT_stop']=tmp['Start']+expand
    tmp_profile=getncRNABindingFrac(tmp)
    plt.subplot(plotDim,plotDim,i)
    bins=np.arange(0,1,0.01)
    hist,bins=np.histogram(tmp_profile['frac'],bins=bins)
    hist=np.array(hist/float(tmp_profile['frac'].shape[0]),dtype=float)
    width=0.7*(bins[1]-bins[0])
    center=(bins[:-1] + bins[1:])/2
    plt.bar(center,hist,align='center',width=width,color='blue',alpha=0.75)
    plt.tick_params(axis='x',labelsize=5) 
    plt.tick_params(axis='y',labelsize=5)  
    plt.xlabel('Fraction of gene body (5p - 3p)',fontsize=5)
    plt.title('Binding profile for %s'%name,fontsize=5)
    i+=1
fig6.tight_layout()
fig6.savefig(outfilepath+'Figure6.png',format='png',bbox_inches='tight',dpi=150,pad_inches=0.5)
fig6.savefig(outfilepath+'Figure6.pdf',format='pdf',bbox_inches='tight',dpi=150,pad_inches=0.5)

# <codecell>

logOpen.close()

# <codecell>



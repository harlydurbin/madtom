#!/bin/python

import sys
import os
import subprocess

SPECIES = sys.argv[1]
ID = sys.argv[2]
LIBRARY = sys.argv[3]		#so far only works for one library, but could change to pass it a comma separated list (AP,BP,CP) then split the list and loop thru cmd1, cmd2, cmd3 for each library
FASTQ_PATH_1 = sys.argv[4]
FASTQ_PATH_2 = sys.argv[5]
REF_GENOME_PATH = sys.argv[6]
PREFIX = SPECIES + '.' + ID + '.' + LIBRARY
LOG_FILE = PREFIX + '.cmdlog'
LOG = open(LOG_FILE, 'w')

###########    ADAPTOR AND QUALITY TRIMMING    ###########
#6 threads optimal; specify phred scale 33; generates paired files plus forward and reverse unpaired (if mate doesn't pass qc)
#drop read if <35 bp after trimming or if avg qual of read <20; trim base if qual <20
cmd1 = 'java -Djava.io.tmpdir=/data/lwwvd/processing/' + ID + '.tmp -XX:ParallelGCThreads=6 -jar /usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 6 -phred33 ' + FASTQ_PATH_1 + ' ' + FASTQ_PATH_2 + ' ' + PREFIX + '.1.P.fastq ' + SPECIES + '.' + ID + '.' + LIBRARY + '.1.U.fastq ' + PREFIX + '.2.P.fastq ' + PREFIX + '.2.U.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 AVGQUAL:20 MINLEN:35 ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.35/adapters/TruSeq3-PE.fa'
LOG.write(cmd1 + "\n")
subprocess.call(cmd1, shell=True)

###########    BWA ALIGNMENTS    ###########
#align paired files (cmd2p), unpaired 1 (cmd2u1), and unpaired 2 (cmd2u2) from trimmomatic
#generate picard compatible files/read group names
#use 62 threads (where available/on mugenomics6)
cmd2p = 'bwa mem -M -t 62 -R \'@RG\\tID:' + ID + '\\tSM:' + SPECIES + ID + '\\tLB:a\\tPL:ILLUMINA\' ' + REF_GENOME_PATH + ' ' + PREFIX + '.1.P.fastq ' + PREFIX + '.2.P.fastq > ' + PREFIX + '.P.sam'
LOG.write(cmd2p + "\n")
subprocess.call(cmd2p, shell=True)

cmd2u1 = 'bwa mem -M -t 62 -R \'@RG\\tID:' + ID + '\\tSM:' + SPECIES + ID + '\\tLB:a\\tPL:ILLUMINA\' ' + REF_GENOME_PATH + ' ' + PREFIX + '.1.U.fastq > ' + PREFIX + '.U1.sam'
LOG.write(cmd2u1 + "\n")
subprocess.call(cmd2u1, shell=True)

cmd2u2 = 'bwa mem -M -t 62 -R \'@RG\\tID:' + ID + '\\tSM:' + SPECIES + ID + '\\tLB:a\\tPL:ILLUMINA\' ' + REF_GENOME_PATH + ' ' + PREFIX + '.2.U.fastq > ' + PREFIX + '.U2.sam'
LOG.write(cmd2u2 + "\n")
subprocess.call(cmd2u2, shell=True)

##########    CREATE MERGED BAM    ##########
#sort paired files (cmd3p), unpaired 1 (cmd3u1), and unpaired 2 (cmd3u2) and output as BAM
cmd3p = 'samtools sort -l 0 -m 100G -\@ 4 -o ' + PREFIX + '.P.sorted.bam -T ' + PREFIX + '.P.sam.TMP ' + PREFIX + '.P.sam'
LOG.write(cmd3p + "\n")
subprocess.call(cmd3p, shell=True)

cmd3u1 = 'samtools sort -l 0 -m 100G -\@ 4 -o ' + PREFIX + '.U1.sorted.bam -T ' + PREFIX + '.U1.sam.TMP ' + PREFIX + '.U1.sam'
LOG.write(cmd3u1 + "\n")
subprocess.call(cmd3u1, shell=True)

cmd3u2 = 'samtools sort -l 0 -m 100G -\@ 4 -o ' + PREFIX + '.U2.sorted.bam -T ' + PREFIX + '.U2.sam.TMP ' + PREFIX + '.U2.sam'
LOG.write(cmd3u2 + "\n")
subprocess.call(cmd3u2, shell=True)

#merge bam files; do not compress
cmd4 = 'java -Djava.io.tmpdir=/data/lwwvd/processing/' + ID + '.tmp -XX:ParallelGCThreads=6 -Xmx10g -jar /usr/local/bin/picard-tools-1.138/picard.jar MergeSamFiles INPUT=' + PREFIX + '.P.sorted.bam INPUT=' + PREFIX + '.U1.sorted.bam INPUT=' + PREFIX + '.U2.sorted.bam OUTPUT=' + PREFIX + '.merged.bam USE_THREADING=TRUE MERGE_SEQUENCE_DICTIONARIES=TRUE COMPRESSION_LEVEL=0'
LOG.write(cmd4 + "\n")
subprocess.call(cmd4, shell=True)

##########    MARK/REMOVE DUPLICATES    ##########
#all records are written to bam file with duplicates flagged; add REMOVE_DUPLICATES=true to not output duplicates with flags
cmd5 = 'java -Djava.io.tmpdir=/data/lwwvd/processing/' + ID + '.tmp -XX:ParallelGCThreads=6 -Xmx20g -jar /usr/local/bin/picard-tools-1.138/picard.jar MarkDuplicates INPUT=' + PREFIX + '.merged.bam OUTPUT=' + PREFIX + '.bam METRICS_FILE=' + PREFIX + '.DUP.METRICS MAX_RECORDS_IN_RAM=1000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 COMPRESSION_LEVEL=0'
LOG.write(cmd5 + "\n")
subprocess.call(cmd5, shell=True)

###########    INDEX BAM    ##########
cmd6 = 'samtools index ' + PREFIX + '.bam'
LOG.write(cmd6 + "\n")
subprocess.call(cmd6, shell=True)

###########    COMPRESS FILES    ##########
gzfastq = 'pigz -p 9 *.fastq'
#subprocess.call(gzfastq, shell=True)

###########    REMOVE UNNECESSARY FILES    ##########
rmtmp = 'rm -r /data/lwwvd/processing/' + ID + '.tmp'		#remove java tmp directory
subprocess.call(rmtmp, shell=True)
rmsam = 'rm *.sam'
subprocess.call(rmsam, shell=True)
rmmerged = 'rm *.sorted.bam'
subprocess.call(rmmerged, shell=True)

###########    RUN FASTQC ON TRIMMED READS     ##########
fastqc = 'fastqc ' + PREFIX + '.1.P.fastq.gz ' + PREFIX + '.2.p.fastq.gz'
LOG.write(fastqc + "\n")
subprocess.call(fastqc, shell=True)

###########    RUN STATISTICS ON BAM FILE     ##########
stats = 'samtools flagstat ' + PREFIX + '.bam'
subprocess.call(stats, shell=True)
stats2 = 'java -Djava.io.tmpdir=/data/lwwvd/tmp -XX:ParallelGCThreads=2 -Xmx10g -jar /usr/local/bin/picard-tools-2.1.0/picard.jar CollectMultipleMetrics INPUT=' + PREFIX + '.bam ASSUME_SORTED=true OUTPUT=' + PREFIX + '.bam.metrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=LIBRARY REFERENCE_SEQUENCE=' + REF_GENOME_PATH + ' PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics STOP_AFTER=30000000'
LOG.write(stats + "\n" + stats2 + "\n")
subprocess.call(stats2, shell=True)

###########    SYNC FINAL BAM FILES TO MUG01_N    ##########
cmdsync = 'rsync -v ' + PREFIX + '.* /CIFS/MUG01_N/taylorj/lwwvd/madtom/BAM'

LOG.close()	#output file with commands used

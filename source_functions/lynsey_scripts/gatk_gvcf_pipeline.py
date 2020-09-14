#!/bin/python

import sys
import os
import subprocess
import time
from multiprocessing.pool import ThreadPool

#Produces GVCF files from BAM alignments

ID = sys.argv[1]
BAM_PATH = sys.argv[2]
REF_GENOME_PATH = sys.argv[3]
IntervalFilePrefix = sys.argv[4]		#chromosomes or unplaced scaffolds intervals
n_start_IntervalFiles = sys.argv[5]
n_end_IntervalFiles = sys.argv[6]

LogFile = ID + '.gatk.cmd.LOG'
TimeFile = ID + '.gatk.cmd.TIME'
LOG = open(LogFile, 'w')
TIME = open(TimeFile, 'w')

n_start = int(n_start_IntervalFiles)
n_end = int(n_end_IntervalFiles) + 1

##########     Realigner Target Creator     ##########
start_time = time.time()
def rtc(rtc_cmd):
  subprocess.call(rtc_cmd, shell=True)
  LOG.write(rtc_cmd + "\n")

tp = ThreadPool(14)
for i in range(n_start, n_end):
  interval = str(i)
  intervalfile = IntervalFilePrefix + '.' + interval + '.interval_list'
  rtc_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/' + ID + '.tmp -XX:ParallelGCThreads=6 -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R ' + REF_GENOME_PATH + ' -L ' + intervalfile + ' -I ' + BAM_PATH + ' -o ' + ID + '.' + interval + '.forIndelRealigner.intervals'
  tp.apply_async(rtc, (rtc_cmd,))
  
tp.close()
tp.join()
end_time = time.time()
time_taken = end_time - start_time
TIME.write('RealignerTargetCreator ' + str(round(time_taken, 2)) + '\n')

#########     Indel Realigner    ##########
start_time = time.time()
def ir(ir_cmd):
  subprocess.call(ir_cmd, shell=True)
  LOG.write(ir_cmd + "\n")

tp = ThreadPool(16)
for i in range(n_start, n_end):
  interval = str(i)
  intervalfile = IntervalFilePrefix + '.' + interval + '.interval_list'
  ir_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/' + ID + '.tmp -XX:ParallelGCThreads=6 -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R ' + REF_GENOME_PATH + ' -I ' + BAM_PATH + ' -T IndelRealigner -targetIntervals ' + ID + '.' + interval + '.forIndelRealigner.intervals -L ' + intervalfile + ' -o ' + ID + '.' + interval + '.realigned.bam --bam_compression 0'
  tp.apply_async(ir, (ir_cmd,))

tp.close()
tp.join()
end_time = time.time()
time_taken = end_time - start_time
TIME.write('IndelRealigner ' + str(round(time_taken, 2)) + '\n')

##########     Merge Realigned Files     ##########
start_time = time.time()
filestomerge = []
for i in range(n_start, n_end):
  filename = ID + '.' + str(i) + '.realigned.bam'
  filestomerge.append(filename)
filelist = ID + '.tmp_MergeRealignedFiles.list'
fh = open(filelist, 'w')
for x in filestomerge:
  fh.write(x + "\n")
fh.close()

mrf = 'samtools merge -@ 10 -u -f -c -b ' + ID + '.tmp_MergeRealignedFiles.list ' + ID + '.realigned.bam'
bai = 'samtools index ' + ID + '.realigned.bam'
subprocess.call(mrf, shell=True)
subprocess.call(bai, shell=True)
LOG.write(mrf + "\n")
LOG.write(bai + "\n")
end_time = time.time()
time_taken = end_time - start_time
TIME.write('MergeRealignedFiles ' + str(round(time_taken, 2)) + '\n')

##########     Depth of Coverage     ##########
start_time = time.time()
def doc(doc_cmd):
  subprocess.call(doc_cmd, shell=True)
  LOG.write(doc_cmd + "\n")

tp = ThreadPool(12)
for i in range(n_start, n_end):
  interval = str(i)
  intervalfile = IntervalFilePrefix + '.' + interval + '.interval_list'
  doc_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/' + ID + '.tmp -XX:ParallelGCThreads=6 -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -nt 8 -R ' + REF_GENOME_PATH + ' -I ' + ID + '.realigned.bam -T DepthOfCoverage -L ' + intervalfile + ' -o ' + ID + '.' + interval + '.realigned.bam.coverage -omitBaseOutput -omitIntervals --omitLocusTable'
  tp.apply_async(doc, (doc_cmd,))

tp.close()
tp.join()
end_time = time.time()
time_taken = end_time - start_time
TIME.write('DepthOfCoverage ' + str(round(time_taken, 2)) + '\n')

##########     Haplotype Caller     ##########
start_time = time.time()
def hapcall(hap_cmd):
  subprocess.call(hap_cmd, shell=True)
  LOG.write(hap_cmd + "\n")

tp = ThreadPool(18)
for i in range(n_start, n_end):
  interval = str(i)
  intervalfile = IntervalFilePrefix + '.' + interval + '.interval_list'
  hap_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/' + ID + '.tmp -XX:ParallelGCThreads=4 -Xmx15g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -nct 6 -ERC GVCF -T HaplotypeCaller -R ' + REF_GENOME_PATH + ' -L ' + intervalfile + ' -I ' + ID + '.realigned.bam -o ' + ID + '.' + interval + '.g.vcf.gz --heterozygosity 0.001 --pcr_indel_model NONE'
  tp.apply_async(hapcall, (hap_cmd,))
  
tp.close()
tp.join()
end_time = time.time()
time_taken = end_time - start_time
TIME.write('HaplotypeCaller ' + str(round(time_taken, 2)) + '\n')

#########     END     ##########

LOG.close()
TIME.close()

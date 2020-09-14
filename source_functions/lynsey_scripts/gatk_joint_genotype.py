#!/bin/python

import sys
import os
import subprocess
import time
from multiprocessing.pool import ThreadPool

IDs_to_combine_file = sys.argv[1]		#text file with one ID per line (ending in g.vcf.list)
REF_GENOME_PATH = sys.argv[2]
INTERVAL_FILE_PREFIX = sys.argv[3]
OUT_PREFIX = sys.argv[4]
LOG_FILE = sys.argv[5]
nIntervalFiles = sys.argv[6]	#the total number of interval files

nInt = int(nIntervalFiles)

IDs_fh = open(IDs_to_combine_file, 'r')
IDs = []
for ind in IDs_fh:
  ind = ind.strip()
  IDs.append(ind)

LOG = open(LOG_FILE, 'w')

##########     Genotype GVCFs     ##########
def geno(geno_cmd):
  subprocess.call(geno_cmd, shell=True)
  LOG.write(geno_cmd + "\n")

tp = ThreadPool(6)
for i in range(0,nInt):
  interval = str(i)
  intervalfile = INTERVAL_FILE_PREFIX+ '.' + interval + '.interval_list'
  gvcf_string = ''
  for x in IDs:
    gvcf_string = gvcf_string + ' -V /CIFS/MUG01_N/taylorjerr/lwwvd/madtom/GVCFs/' + x + '/' + x + '.' + interval + '.g.vcf.gz'
  geno_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/genoGVCF.' + interval + '.tmp -XX:ParallelGCThreads=4 -Xmx50g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ' + REF_GENOME_PATH + ' -L ' + intervalfile + ' ' + gvcf_string + ' -o ' + OUT_PREFIX + '.combined.' + interval + '.vcf.gz'
  tp.apply_async(geno, (geno_cmd,))

tp.close()
tp.join()

IDs_fh.close()
LOG.close()

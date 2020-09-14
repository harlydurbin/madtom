#!/bin/python

import sys
import subprocess
from multiprocessing.pool import ThreadPool

REF_GENOME_PATH = sys.argv[1]
JOINT_VCF_FILE_PREFIX = sys.argv[2]
LOG_FILE = sys.argv[3]
N_INTERVALS = sys.argv[4]		#number of intervals, madtom=29 (0 thru 28)
LOG = open(LOG_FILE, 'w')
nInt = int(N_INTERVALS)

##########     SELECT BIALLELIC SNPS     ##########

def snp(snp_cmd):
  subprocess.call(snp_cmd, shell=True)
  LOG.write(snp_cmd + "\n")

tp = ThreadPool(12)
for i in range(0,nInt):
  interval = str(i)
  vcf_file = JOINT_VCF_FILE_PREFIX + '.' + interval + '.vcf.gz'
  snp_vcf_out = JOINT_VCF_FILE_PREFIX + '.' + interval + '.SNP.vcf.gz'
  snp_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/selectSNP.' + interval + '.tmp -XX:ParallelGCThreads=4 -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R ' + REF_GENOME_PATH + ' -V ' + vcf_file + ' -o ' + snp_vcf_out + ' -selectType SNP --restrictAllelesTo BIALLELIC'
  tp.apply_async(snp, (snp_cmd,))

tp.close()
tp.join()

##########     SELECT INDELS     ##########

def indel(indel_cmd):
  subprocess.call(indel_cmd, shell=True)
  LOG.write(indel_cmd + "\n")

tp = ThreadPool(12)
for i in range(0,nInt):
  interval = str(i)
  vcf_file = JOINT_VCF_FILE_PREFIX + '.' + interval + '.vcf.gz'
  indel_vcf_out = JOINT_VCF_FILE_PREFIX + '.' + interval + '.INDEL.vcf.gz'
  indel_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/selectINDEL.' + interval + '.tmp -XX:ParallelGCThreads=4 -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R ' + REF_GENOME_PATH + ' -V ' + vcf_file + ' -o ' + indel_vcf_out + ' -selectType INDEL'
  tp.apply_async(indel, (indel_cmd,))

tp.close()
tp.join()

##########     FILTER SNPS     #########

##Filter SNPs on: 
#QD ("Variant Confidence/Quality by Depth"), 
#FS ("Phred-scaled p-value using Fisher's exact test to detect strand bias"), 
#SOR ("Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"), 
#MQ ("RMS Mapping Quality"), 
#MQRankSum ("Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"), 
#ReadPosRankSum ("Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias")

def fs(fs_cmd):
  subprocess.call(fs_cmd, shell=True)
  LOG.write(fs_cmd + "\n")

tp = ThreadPool(12)
for i in range(0,nInt):
  interval = str(i)
  snp_vcf = JOINT_VCF_FILE_PREFIX + '.' + interval + '.SNP.vcf.gz'
  filtered_output = JOINT_VCF_FILE_PREFIX + '.' + interval + '.filter.SNP.vcf.gz'
  fs_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/filterVCF_SNP.' + interval + '.tmp -XX:ParallelGCThreads=6 -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF_GENOME_PATH + ' -V ' + snp_vcf + ' -filter "QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -filterName "SNP_filter" -o ' + filtered_output
  tp.apply_async(fs, (fs_cmd,))

tp.close()
tp.join()

##########     FILTER INDELS     ##########

##Filter INDELs on:
#QD, FS, SOR, InbreedingCoeff ("Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"), ReadPosRankSum

def fi(fi_cmd):
  subprocess.call(fi_cmd, shell=True)
  LOG.write(fi_cmd + "\n")

tp = ThreadPool(12)
for i in range(0,nInt):
  interval = str(i)
  indel_vcf = JOINT_VCF_FILE_PREFIX + '.' + interval + '.INDEL.vcf.gz'
  filtered_output = JOINT_VCF_FILE_PREFIX + '.' + interval + '.filter.INDEL.vcf.gz'
  fi_cmd = 'java -Djava.io.tmpdir=/data/lwwvd/filterVCF_INDEL.' + interval + '.tmp -XX:ParallelGCThreads=6 -Xmx24g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T VariantFiltration -R ' + REF_GENOME_PATH + ' -V ' + indel_vcf + ' -filter "QD < 2.0 || FS>200.0 || SOR>10.0 || InbreedingCoeff<-0.8 || ReadPosRankSum < -20.0" -filterName "INDEL_filter" -o ' + filtered_output
  tp.apply_async(fi, (fi_cmd,))

tp.close()
tp.join()

##########     CONCATENATE VCFs     ##########
'''
def cat(cat_vcfs):
  subprocess.call(cat_vcfs, shell=True)
  LOG.write(cat_vcfs)

cat_vcfs = 'java -Djava.io.tmpdir=/data/lwwvd/catVCFs.tmp -cp /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ' + REF_GENOME_PATH + ' -out output.vcf -assumeSorted'
for i in range(3,50):
  interval = str(i)
  cat_vcfs = cat_vcfs + ' -V ' + 
cat(cat_vcfs)
'''

#Can also concatenate (for GWAS): bcftools concat -f filterSNPvcfs.list -O v -o BBUB.combined.filter.SNP.vcf
#then, convert to plink format (also use make_fake_plink_posns.py to correct for too many contig IDs) vcftools --vcf BBUB.combined.filter.SNP.vcf --out BBUB.combined.filter.plink --plink

LOG.close()

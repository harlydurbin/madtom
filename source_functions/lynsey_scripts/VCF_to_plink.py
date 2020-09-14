#!/usr/bin/python

import sys
import subprocess

VCF_FILE_PREFIX = sys.argv[1]
N_INTERVALS = sys.argv[2]		#number of intervals, madtom=29 (0 thru 28)
nInt = int(N_INTERVALS)


##########    CONVERT VCF TO PLINK    ##########

for i in range(0,nInt):
  interval = str(i)
  #plink_cmd = 'plink --vcf ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP.vcf.gz --allow-extra-chr --allow-no-sex --geno 0.10 --make-bed --out ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP' 
  #subprocess.call(plink_cmd, shell=True)

##########     FILTER TO VARIABLE SNPs    #########
#must have minor allele frequence > 1/22

for i in range(0,nInt):
  interval = str(i)
  #plink_cmd2 = 'plink --bfile ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP --allow-extra-chr --allow-no-sex --maf 0.04 --geno 0.10 --make-bed --out ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP.mdtmvar'
  #subprocess.call(plink_cmd2, shell=True)
  plink_cmd3 = 'plink --bfile ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP.mdtmvar --allow-extra-chr --allow-no-sex --remove slender_madtom_id.txt --maf 0.05 --make-bed --out ' + VCF_FILE_PREFIX + '.' + interval + '.filter.SNP.neoshovar'
  subprocess.call(plink_cmd3, shell=True)

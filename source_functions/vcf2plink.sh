#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH -p General,Lewis,BioCompute,htc4,hpc5  # for jobs < 2hrs try 'General'
#SBATCH -c 12 # cores per task
#SBATCH --mem=40G  # memory per core (default is 1GB/core)
#SBATCH --time 01:00:00  # days-hours:minutes
#SBATCH --account=animalsci  # investors will replace this with their account name
#SBATCH --output=log/slurm_out/vcf2plink/vcf2plink.%j.out # %j is the unique jobID

module load plink/plink-1.90b

plink --vcf output/joint_genotyping/madtom.sorted_snps.vcf.gz --allow-extra-chr --geno 0.01 --maf 0.01 --make-bed --threads 6 --out output/faststructure/madtom.sorted_snps

plink --bfile output/faststructure/madtom.sorted_snps --allow-extra-chr --remove output/faststructure/remove_stonecat.txt --make-bed --threads 6 --out output/faststructure/madtom.sorted_snps.no_stonecat

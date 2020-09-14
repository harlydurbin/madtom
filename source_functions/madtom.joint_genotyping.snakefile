# snakemake -s source_functions/madtom.joint_genotyping.snakefile -j 400 --rerun-incomplete --latency-wait 30 --config --cluster-config source_functions/cluster_config/madtom.joint_genotyping.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/200317.madtom.joint_genotyping.log

import pandas as pd
import os

configfile: "source_functions/config/madtom.joint_genotyping.config.yaml"

# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/joint_genotyping", exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

# Make temp directories if they don't exist
for x in expand("temp/joint_genotyping/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

# Import list of intervals, convert to a list in order to expand over
intervals_df = pd.read_csv("output/call_genotypes/interval_list/all_intervals.csv")

all_intervals = intervals_df["group_number"].tolist()

#all_intervals = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

rule all:
	input:
		# expand("output/joint_genotyping/{sample_number}/{sample_number}.{interval}.g.vcf.gz", sample_number = config['sample_number'], interval = all_intervals)
		"output/joint_genotyping/madtom.sorted_snps.vcf.gz", "output/joint_genotyping/madtom.sorted_snps.vcf.gz.tbi"

rule genotype_gvcfs:
	input:
		gvcfs = lambda wildcards: expand("output/call_genotypes/{sample_number}/{sample_number}.{interval}.g.vcf.gz", sample_number = config['sample_number'], interval = wildcards.interval),
		interval_list = "output/call_genotypes/interval_list/{interval}.interval_list"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['genotype_gvcfs_gc'],
		xmx = config['genotype_gvcfs_xmx'],
		java_tmp = "temp/joint_genotyping/genotype_gvcfs/{interval}",
		gatk_path = config['gatk_path'],
		arguments = lambda wildcards: expand("-V output/call_genotypes/{sample_number}/{sample_number}.{interval}.g.vcf.gz", sample_number = config['sample_number'], interval = wildcards.interval),
		psrecord = "log/psrecord/joint_genotyping/genotype_gvcfs.{interval}.log"
	output:
		vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.vcf.gz"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T GenotypeGVCFs -R {params.ref_genome} -L {input.interval_list} {params.arguments} -o {output.vcf}" --log {params.psrecord} --include-children --interval 2
		"""

# Restrict to biallelic SNPs
rule select_snps:
	input:
		vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.vcf.gz"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['select_variants_gc'],
		xmx = config['select_variants_xmx'],
		java_tmp = "temp/joint_genotyping/select_snps/{interval}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/select_snps.{interval}.log"
	output:
		snp_vcf = temp("output/joint_genotyping/{interval}/MDTM.combined.{interval}.SNP.vcf.gz")
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T SelectVariants -R {params.ref_genome} -V {input.vcf} -o {output.snp_vcf} -selectType SNP --restrictAllelesTo BIALLELIC" --log {params.psrecord} --include-children --interval 2
		"""

rule filter_snps:
	input:
		snp_vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.SNP.vcf.gz"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['variant_filtration_gc'],
		xmx = config['variant_filtration_xmx'],
		java_tmp = "temp/joint_genotyping/filter_snps/{interval}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/filter_snps.{interval}.log"
	output:
		filtered_snp_vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.filter.SNP.vcf.gz"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T VariantFiltration -R {params.ref_genome} -V {input.snp_vcf} -filter 'QD < 2.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' -filterName 'SNP_filter' -o {output.filtered_snp_vcf}" --log {params.psrecord} --include-children --interval 2
		"""

# Generate list of paths to VCFs for merge_vcfs, write to a file
# This was making Snakemake take way too long to build a DAG so i just did it in R
rule merge_list:
	input:
		vcfs = expand("output/joint_genotyping/{interval}/MDTM.combined.{interval}.filter.SNP.vcf.gz", interval = all_intervals)
	output:
		merge_list = "output/joint_genotyping/merge_files.list"
	params:
		path_list = lambda wildcards: expand("output/joint_genotyping/{interval}/MDTM.combined.{interval}.filter.SNP.vcf.gz\n", interval = all_intervals)
	run:
		with open(output.merge_list, "w+") as outfile:
			outfile.write(params.path_list)

rule concat_sort:
	input:
		vcfs = expand("output/joint_genotyping/{interval}/MDTM.combined.{interval}.filter.SNP.vcf.gz", interval = all_intervals),
		merge_list = "output/joint_genotyping/merge_files.list"
	params:
		psrecord = "log/psrecord/joint_genotyping/concat_sort.log",
		bcftools_module = config['bcftools_module'],
		nt = config['concat_sort_nt'],
		temp = "temp/joint_genotyping/concat_sort"
	output:
		sorted_vcf = "output/joint_genotyping/madtom.sorted_snps.vcf.gz"
	# -a: allow overlaps: First coordinate of the next file can precede last record of the current file.
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools concat -a -f {input.merge_list} -Ou --threads {params.nt} -f {input.merge_list} | bcftools sort --output-type z --output-file {output.sorted_vcf} --temp-dir {params.temp}" --log {params.psrecord} --include-children --interval 2
		"""

rule index:
	input:
		sorted_vcf = "output/joint_genotyping/madtom.sorted_snps.vcf.gz"
	params:
		psrecord = "log/psrecord/joint_genotyping/index.log",
		samtools_module = config['samtools_module'],
		htslib_module = config['htslib_module']
	output:
		tbi = "output/joint_genotyping/madtom.sorted_snps.vcf.gz.tbi"
	shell:
		"""
		module load {params.htslib_module}
		module load {params.samtools_module}
		psrecord "tabix {input.sorted_vcf}" --log {params.psrecord} --include-children --interval 2
		"""

# Plink logs show Lynsey filtered for MAF and variants with many missing calls
#   --geno 0.01
#  --maf 0.01

#### SCRATCH ####

# 200312 Jared said to exclude indels

# rule select_indels:
# 	input:
# 		vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.vcf.gz"
# 	params:
# 		java_module = config['java_module'],
# 		gatk_module = config['gatk_module'],
# 		ref_genome = config['ref_genome'],
# 		gc_threads = config['select_variants_gc'],
# 		xmx = config['select_variants_xmx'],
# 		java_tmp = "temp/joint_genotyping/select_indels/{interval}",
# 		gatk_path = config['gatk_path'],
# 		psrecord = "log/psrecord/joint_genotyping/select_indels.{interval}.log"
# 	output:
# 		indel_vcf = temp("output/joint_genotyping/{interval}/MDTM.combined.{interval}.INDEL.vcf.gz")
# 	shell:
# 		"""
# 		module load {params.java_module}
# 		module load {params.gatk_module}
# 		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads= -Xmx{params.xmx}g -jar /usr/local/bin/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R {params.ref_genome} -V {input.vcf} -o {output.indel_vcf} -selectType INDEL" --log {params.psrecord} --include-children --interval 2
# 		"""

# rule filter_indels:
# 	input:
# 		indel_vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.INDEL.vcf.gz"
# 	params:
# 		java_module = config['java_module'],
# 		gatk_module = config['gatk_module'],
# 		ref_genome = config['ref_genome'],
# 		gc_threads = config['variant_filtration_gc'],
# 		xmx = config['variant_filtration_xmx'],
# 		java_tmp = "temp/joint_genotyping/filter_indels/{interval}",
# 		gatk_path = config['gatk_path'],
# 		psrecord = "log/psrecord/joint_genotyping/filter_indels.{interval}.log"
# 	output:
# 		filtered_indel_vcf = "output/joint_genotyping/{interval}/MDTM.combined.{interval}.filter.INDEL.vcf.gz"
# 	shell:
# 		"""
# 		module load {params.java_module}
# 		module load {params.gatk_module}
# 		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T VariantFiltration -R {params.ref_genome} -V {input.indel_vcf} -filter "QD < 2.0 || FS>200.0 || SOR>10.0 || InbreedingCoeff<-0.8 || ReadPosRankSum < -20.0" -filterName 'INDEL_filter' -o {output.filtered_indel_vcf}" --log {params.psrecord} --include-children --interval 2
# 		"""

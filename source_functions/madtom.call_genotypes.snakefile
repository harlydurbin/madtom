# snakemake -s source_functions/madtom.call_genotypes.snakefile -j 800 --rerun-incomplete --latency-wait 60 --config --cluster-config source_functions/cluster_config/madtom.call_genotypes.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account}" -p &> log/snakemake_log/200316.madtom.call_genotypes.log

import pandas as pd
import os

configfile: "source_functions/config/madtom.call_genotypes.config.yaml"

# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

# Import list of intervals, convert to a list in order to expand over
intervals_df = pd.read_csv("output/call_genotypes/interval_list/all_intervals.csv")

all_intervals = intervals_df["group_number"].tolist()

#all_intervals = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

rule all:
	input:
		expand("output/call_genotypes/{sample_number}/{sample_number}.{interval}.g.vcf.gz", sample_number = config['sample_number'], interval = all_intervals)

rule create_dict:
	input:
		ref_genome = config['ref_genome']
	params:
		java_module = config['java_module'],
		picard_module = config['picard_module'],
		picard_path = config['picard_path']
	output:
		dict = config['dict']
	shell:
		"""
		module load {params.java_module}
		module load {params.picard_module}
		java -jar {params.picard_path} CreateSequenceDictionary R= {input.ref_genome} O= {output.dict}
		"""

rule create_fai:
	input:
		ref_genome = config['ref_genome']
	params:
		samtools_module = config['samtools_module']
	output:
		fai = config['ref_genome']+'.fai'
	shell:
		"""
		module load {params.samtools_module}
		samtools faidx {input.ref_genome}
		"""

rule target_creator:
	input:
		dict = config['dict'],
		fai = config['ref_genome']+'.fai',
		interval_list = "output/call_genotypes/interval_list/{interval}.interval_list",
		#indexed_bam = "avail.txt"
		indexed_bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam",
		bai = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam.bai"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['target_creator_gc'],
		xmx = config['target_creator_xmx'],
		java_tmp = "temp/call_genotypes/{sample_number}",
		gatk_path = config['gatk_path'],
		threads = config['target_creator_threads'],
		psrecord = "log/psrecord/call_genotypes/target_creator.{sample_number}.{interval}.log"
	output:
		intervals = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.forIndelRealigner.intervals"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -nt {params.threads} -T RealignerTargetCreator -R {params.ref_genome} -L {input.interval_list} -I {input.indexed_bam} -o {output.intervals}" --log {params.psrecord} --include-children --interval 2
		"""

rule indel_realigner:
	input:
		indexed_bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam",
		bai = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam.bai",
		interval_list = "output/call_genotypes/interval_list/{interval}.interval_list",
		intervals = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.forIndelRealigner.intervals"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['indel_realigner_gc'],
		xmx = config['indel_realigner_xmx'],
		threads = config['indel_realigner_threads'],
		java_tmp = "temp/call_genotypes/{sample_number}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/call_genotypes/indel_realigner.{sample_number}.{interval}.log"
	output:
		realigned_bam = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -R {params.ref_genome} -I {input.indexed_bam} -T IndelRealigner -targetIntervals {input.intervals} -L {input.interval_list} -o {output.realigned_bam} --bam_compression 0" --log {params.psrecord} --include-children --interval 2
		"""
rule merge_realigned:
	input:
		realigned_bam = expand("output/call_genotypes/{{sample_number}}/{{sample_number}}.{interval}.realigned.bam", interval = all_intervals)
	params:
		realigned_wildcard = "output/call_genotypes/{sample_number}/{sample_number}.*.realigned.bam",
		threads = config['indel_merge_threads'],
		psrecord = "log/psrecord/call_genotypes/merge_realigned.{sample_number}.log",
		samtools_module = config['samtools_module']
	output:
		merge_list = "output/call_genotypes/{sample_number}/{sample_number}.tmp_MergeRealignedFiles.list",
		merged_bam = "output/call_genotypes/{sample_number}/{sample_number}.realigned.bam"
	shell:
		"""
		module load {params.samtools_module}
		psrecord "ls {params.realigned_wildcard} | tr '\t' '\n' > {output.merge_list};
		samtools merge -@ {params.threads} -u -f -c -b {output.merge_list} {output.merged_bam};
		samtools index {output.merged_bam}" --log {params.psrecord} --include-children --interval 2
		"""

rule depth_of_coverage:
	input:
		merged_bam = "output/call_genotypes/{sample_number}/{sample_number}.realigned.bam",
		interval_list = "output/call_genotypes/interval_list/{interval}.interval_list"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['depth_of_coverage_gc'],
		xmx = config['depth_of_coverage_xmx'],
		threads = config['depth_of_coverage_threads'],
		java_tmp = "temp/call_genotypes/{sample_number}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/call_genotypes/depth_of_coverage.{sample_number}.{interval}.log",
		# output prefix
		bam_coverage = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam.coverage"
	output:
		stat = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam.coverage.sample_statistics",
		sum = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam.coverage.sample_summary"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -nt {params.threads} -R {params.ref_genome} -I {input.merged_bam} -T DepthOfCoverage -L {input.interval_list} -o {params.bam_coverage} -omitBaseOutput -omitIntervals --omitLocusTable" --log {params.psrecord} --include-children --interval 2
		"""

rule haplotype_caller:
	input:
		merged_bam = "output/call_genotypes/{sample_number}/{sample_number}.realigned.bam",
		stat = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam.coverage.sample_statistics",
		sum = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.realigned.bam.coverage.sample_summary",
		interval_list = "output/call_genotypes/interval_list/{interval}.interval_list"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['haplotype_caller_gc'],
		xmx = config['haplotype_caller_xmx'],
		threads = config['haplotype_caller_threads'],
		java_tmp = "temp/call_genotypes/{sample_number}/",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/call_genotypes/haplotype_caller.{sample_number}.{interval}.log"
	output:
		gvcf = "output/call_genotypes/{sample_number}/{sample_number}.{interval}.g.vcf.gz"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -nct {params.threads} -ERC GVCF -T HaplotypeCaller -R {params.ref_genome} -L {input.interval_list} -I {input.merged_bam} -o {output.gvcf} --heterozygosity 0.001 --pcr_indel_model NONE" --log {params.psrecord} --include-children --interval 2
		"""

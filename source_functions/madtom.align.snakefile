# snakemake -s source_functions/madtom.align.snakefile -j 80 --rerun-incomplete --latency-wait 60 --config --cluster-config source_functions/cluster_config/madtom.align.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/200306.madtom.align.log

configfile: "source_functions/config/madtom.align.config.yaml"

rule target:
	input:
		target = expand("output/align/{sample_number}/MDTM.{sample_number}.aa.bam.bai", sample_number = config['sample_number'])

rule index_reference:
	input:
		ref_genome = config['ref_genome']
	params:
		bwa_module = config['bwa_module'],
		psrecord = "log/psrecord/align/index_reference.log"
	output:
		fm_index = expand("{ref_genome}.bwt", ref_genome = config['ref_genome'])
	shell:
		"""
		module load {params.bwa_module}
		psrecord "bwa index {input.ref_genome}" --log {params.psrecord} --include-children --interval 2
		"""

# Input is output of trimmomatic
# once for paired {1, 2}
rule bwa_mem_paired:
	input:
		fastq = expand("/storage/htc/deckerlab/FROM_MUG01/MUG01_N/lwwvd/madtom/BAM/ipcoco_ref/{{sample_number}}/MDTM.{{sample_number}}.aa.{read_num}.P.fastq.gz", read_num = config['read_num']),
		fm_index = expand("{ref_genome}.bwt", ref_genome = config['ref_genome'])
	params:
		# Denovo assembly
		ref_genome = config['ref_genome'],
		threads = config['bwa_threads'],
		bwa_module = config['bwa_module'],
		# Read group header
		rg_header = "'@RG\\tID:{sample_number}\\tSM:MDTM{sample_number}\\tLB:a\\tPL:ILLUMINA'",
		psrecord = "log/psrecord/align/bwa_mem_paired.{sample_number}.log"
	output:
		sam = "output/align/{sample_number}/MDTM.{sample_number}.aa.P.sam",
	shell:
		"""
		module load {params.bwa_module}
		psrecord "bwa mem -M -t {params.threads} -R {params.rg_header} {params.ref_genome} {input.fastq} > {output.sam}" --log {params.psrecord} --include-children --interval 2
		"""

# once for unpaired {1}
# once for unpaired {2}
rule bwa_mem_unpaired:
	input:
		fastq = "/storage/htc/deckerlab/FROM_MUG01/MUG01_N/lwwvd/madtom/BAM/ipcoco_ref/{sample_number}/MDTM.{sample_number}.aa.{read_num}.U.fastq.gz",
		fm_index = expand("{ref_genome}.bwt", ref_genome = config['ref_genome'])
	params:
		# Denovo assembly
		ref_genome = config['ref_genome'],
		threads = config['bwa_threads'],
		bwa_module = config['bwa_module'],
		# Read group header
		rg_header = "'@RG\\tID:{sample_number}\\tSM:MDTM{sample_number}\\tLB:a\\tPL:ILLUMINA'",
		psrecord = "log/psrecord/align/bwa_mem_unpaired.{sample_number}.U{read_num}.log"
	output:
		sam = "output/align/{sample_number}/MDTM.{sample_number}.aa.U{read_num}.sam"
	shell:
		"""
		module load {params.bwa_module}
		psrecord "bwa mem -M -t {params.threads} -R {params.rg_header} {params.ref_genome} {input.fastq} > {output.sam}" --log {params.psrecord} --include-children --interval 2
		"""

# samtools sort -l 0 -m 100G -\@ 4 -o MDTM.93678.aa.U2.sorted.bam -T MDTM.93678.aa.U2.sam.TMP MDTM.93678.aa.U2.sam
rule samtools_sort:
	input:
		sam = "output/align/{sample_number}/MDTM.{sample_number}.aa.{read_id}.sam"
	params:
		samtools_module = config['samtools_module'],
		samtools_mem = config['samtools_mem'],
		samtools_threads = config['samtools_threads'],
		sam_tmp = "output/align/{sample_number}/MDTM.{sample_number}.aa.{read_id}.sam.TMP",
		psrecord = "log/psrecord/align/samtools_sort.{sample_number}.{read_id}.log"
	output:
		bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.{read_id}.sorted.bam"
	shell:
		"""
		module load {params.samtools_module}
		psrecord "samtools sort -l 0 -m {params.samtools_mem}G -@ {params.samtools_threads} -o {output.bam} -T {params.sam_tmp} {input.sam}"  --log {params.psrecord} --include-children --interval 2
		"""

rule merge_bam:
	input:
		bam = expand("output/align/{{sample_number}}/MDTM.{{sample_number}}.aa.{read_id}.sorted.bam", read_id = config['read_id'])
	params:
		java_tmp = "temp/align/{sample_number}",
		gc_threads = config['picard_gc'],
		xmx = config['merge_bam_xmx'],
		picard_path = config['picard_path'],
		picard_module = config['picard_module'],
		java_module = config['java_module'],
		input_statement = lambda wildcards: expand("INPUT=output/align/{sample_number}/MDTM.{sample_number}.aa.{read_id}.sorted.bam", sample_number = wildcards.sample_number, read_id = config['read_id']),
		psrecord = "log/psrecord/align/merge_bam.{sample_number}.log"
	output:
		merged_bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.merged.bam"
	shell:
		"""
		module load {params.java_module}
		module load {params.picard_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.picard_path} MergeSamFiles {params.input_statement} OUTPUT={output.merged_bam} USE_THREADING=TRUE MERGE_SEQUENCE_DICTIONARIES=TRUE COMPRESSION_LEVEL=0" --log {params.psrecord} --include-children --interval 2
		"""

rule mark_duplicates:
	input:
		merged_bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.merged.bam"
	params:
		java_tmp = "temp/align/{sample_number}",
		gc_threads = config['picard_gc'],
		xmx = config['mark_duplicates_xmx'],
		picard_path = config['picard_path'],
		picard_module = config['picard_module'],
		java_module = config['java_module'],
		metrics_file = "output/align/{sample_number}/MDTM.{sample_number}.aa.DUP.METRICS",
		psrecord = "log/psrecord/align/mark_duplicates.{sample_number}.log"
	output:
		bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam"
	shell:
		"""
		module load {params.java_module}
		module load {params.picard_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.picard_path} MarkDuplicates INPUT={input.merged_bam} OUTPUT={output.indexed_bam} METRICS_FILE={params.metrics_file} MAX_RECORDS_IN_RAM=1000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 COMPRESSION_LEVEL=0" --log {params.psrecord1} --include-children --interval 20 --log {params.psrecord}
		psrecord "samtools index {output.indexed_bam}" --log {params.psrecord2} --include-children --interval 2
		"""

rule index_bam:
	input:
		bam = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam"
	params:
		samtools_module = config['samtools_module'],
		psrecord = "log/psrecord/align/index_bam.{sample_number}.log"
	output:
		bai = "output/align/{sample_number}/MDTM.{sample_number}.aa.bam.bai"
	shell:
		"""
		module load {params.samtools_module}
		psrecord "samtools index {input.bam}" --log {params.psrecord} --include-children --interval 2
		"""

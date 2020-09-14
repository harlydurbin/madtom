# snakemake -s source_functions/madtom.faststructure.snakefile -j 4 --rerun-incomplete --latency-wait 30 --config -p &> log/snakemake_log/200325.madtom.faststructure.log

configfile: "source_functions/config/madtom.faststructure.config.yaml"

fsdir = "output/faststructure/"

# Make log directories if they don't exist
# for x in expand("log/slurm_out/{rules}", rules = config['rules']):
#     os.makedirs(x, exist_ok = True)
#
# os.makedirs("log/psrecord/faststructure", exist_ok = True)

rule all:
	input:
		expand("{fsdir}{prefix}.ML.txt", prefix = config['prefix'], fsdir = fsdir)


rule structure:
	input:
		bed = fsdir + "{prefix}.bed",
		bim = fsdir + "{prefix}.bim",
		fam = fsdir + "{prefix}.fam"
	params:
		prefix = fsdir + "{prefix}",
		pop = "{k}",
		faststructure_path = config['faststructure_path']
		# psrecord = "log/psrecord/faststructure/structure.{prefix}.{k}.log"
	output:
		meandP = fsdir + "{prefix}.{k}.meanP",
		meandQ = fsdir + "{prefix}.{k}.meanQ",
		varP = fsdir + "{prefix}.{k}.varP",
		varQ = fsdir + "{prefix}.{k}.varQ"
	shell:
		"""
		/usr/local/bin/anaconda/bin/python {params.faststructure_path}/structure.py -K {params.pop} --input={params.prefix} --output={params.prefix} --prior=simple --cv=0 --full
		"""


rule structure_ml:
	input:
		meandP =
		lambda wildcards:  expand("{fsdir}{prefix}.{k}.meanP", prefix = wildcards.prefix, k = config['k'], fsdir = fsdir),
		meandQ =
		lambda wildcards:  expand("{fsdir}{prefix}.{k}.meanQ", prefix = wildcards.prefix, k = config['k'], fsdir = fsdir),
		varP =
		lambda wildcards:  expand("{fsdir}{prefix}.{k}.varP", prefix = wildcards.prefix, k = config['k'], fsdir = fsdir),
		varQ =
		lambda wildcards:  expand("{fsdir}{prefix}.{k}.varQ", prefix = wildcards.prefix, k = config['k'], fsdir = fsdir),
	params:
		faststructure_path = config['faststructure_path'],
		prefix = fsdir + "{prefix}"
	output:
		ML = fsdir + "{prefix}.ML.txt"

	shell:
		"""
		/usr/local/bin/anaconda/bin/python {params.faststructure_path}/chooseK.py --input={params.prefix} &> {output.ML}
		"""

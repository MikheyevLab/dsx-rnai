# dsx rnai
from snakemake.utils import R

outDir = "/work/MikheyevU/sasha/dsx-rnai/data"   # symlink to /work
refDir = "ref" 	  # symlink to /work
SCRATCH  = "/work/scratch/sasha"

kallistoIndex = "ref/hbee"


GALBRAITH = {"SRR3033273": "active", "SRR3033272": "active", "SRR3033271": "active", "SRR3033270": "active", "SRR3033269": "active", "SRR3033268": "active", "SRR3033267": "active", "SRR3033266": "active", "SRR3033265": "active", "SRR3033264": "active", "SRR3033263": "active", "SRR3033262": "active", "SRR3033261": "inactive", "SRR3033260": "inactive", "SRR3033259": "inactive", "SRR3033258": "inactive", "SRR3033257": "inactive", "SRR3033256": "inactive", "SRR3033255": "inactive", "SRR3033254": "inactive", "SRR3033253": "inactive", "SRR3033252": "inactive", "SRR3033251": "inactive", "SRR3033250": "inactive", "SRR3033249": "inactive"}
HE = {"SRR3102934" : "2-day old worker larvae", "SRR3123272" : "2-day old worker larvae", "SRR3123273" : "2-day old worker larvae", "SRR3123275" : "2-day old worker larvae", "SRR3123276" : "2-day old worker larvae", "SRR3123277" : "2-day old worker larvae", "SRR3123279" : "4-day old worker larvae", "SRR3123281" : "4-day old worker larvae", "SRR3123337" : "4-day old worker larvae", "SRR3123340" : "4-day old worker larvae", "SRR3123341" : "4-day old worker larvae", "SRR3123342" : "4-day old worker larvae", "SRR3123355" : "2-day old queen larvae", "SRR3123357" : "2-day old queen larvae", "SRR3123359" : "2-day old queen larvae", "SRR3123361" : "2-day old queen larvae", "SRR3123362" : "2-day old queen larvae", "SRR3123364" : "2-day old queen larvae", "SRR3123372" : "4-day old queen larvae", "SRR3123380" : "4-day old queen larvae", "SRR3123385" : "4-day old queen larvae", "SRR3123388" : "4-day old queen larvae", "SRR3123389" : "4-day old queen larvae", "SRR3123390" : "4-day old queen larvae", "SRR3123400" : "2-day old drone larvae", "SRR3123402" : "2-day old drone larvae", "SRR3123404" : "2-day old drone larvae", "SRR3123406" : "2-day old drone larvae", "SRR3123407" : "2-day old drone larvae", "SRR3123408" : "2-day old drone larvae", "SRR3123443" : "4-day old drone larvae", "SRR3123445" : "4-day old drone larvae", "SRR3123446" : "4-day old drone larvae", "SRR3123448" : "4-day old drone larvae", "SRR3123449" : "4-day old drone larvae", "SRR3123451" : "4-day old drone larvae"}

SAMPLES = [11, 31, 36, 42, 46, 4, 54, 55, 61, 7, 10, 13, 14, 1, 25, 26, 45, 53, 57, 9] # first 10 GFP
localrules: all, collectRsem, rsemEbseq, collectRsemGalbraith, collectRsemHe

rule all:
	input: outDir + "/rsem/genes.ebseq", expand(outDir + "/kallisto/{sample}/abundance.tsv", sample = SAMPLES)

rule getGalbraith:
	output: outDir + "/galbraith/raw/{accession}_1.fastq.gz", outDir + "/galbraith/raw/{accession}_2.fastq.gz"
	resources: mem=5, time=60*25*3
	shell:
		"""
		module load sra-tools/2.8.2-1
		fastq-dump --gzip --split-files  --outdir {outDir}/galbraith/raw {wildcards.accession}
		"""

rule trimGalbraith:
	input: 
		read1 = outDir + "/galbraith/raw/{accession}_1.fastq.gz",
		read2 = outDir + "/galbraith/raw/{accession}_2.fastq.gz"
	output:
		read1 = outDir + "/galbraith/reads_trimmed/{accession}_R1.fastq.gz",
		read2 = outDir + "/galbraith/reads_trimmed/{accession}_R2.fastq.gz",
	resources: mem=10, time=60*24
	threads: 6
	shell:
		"""
		module load adapterremoval/2.2.2
		AdapterRemoval --trimwindows 5 --minquality 20 --file1 {input.read1} --file2 {input.read2} --output1 {output.read1} --output2 {output.read2} --basename {outDir}/galbraith/reads_trimmed/{wildcards.accession} --gzip --threads {threads}
		"""

rule rsemGalbraith:
	input: rules.trimGalbraith.output
	output: outDir + "/galbraith/rsem/{accession}.genes.results"
	resources: mem=80, time=60*4
	threads: 12
	shell:
		"""
		module load rsem bowtie
		cd {outDir}/galbraith/rsem
		rsem-calculate-expression -p {threads} --paired-end <(zcat {input[0]}) <(zcat {input[1]}) {outDir}/../ref/hbee {wildcards.accession}
		"""

rule collectRsemGalbraith:
	input: expand(outDir + "/galbraith/rsem/{accession}.genes.results", accession = GALBRAITH.keys())
	output: 
		counts = outDir + "/galbraith/rsem/galbraith_genes_counts.csv", 
		tpm = outDir + "/galbraith/rsem/galbraith_tpm.csv"
	shell:
		"""
		python scripts/collect_counts.py genes {outDir}/galbraith/rsem/ > {output.counts}
		python scripts/collect_tpm.py genes {outDir}/galbraith/rsem/ > {output.tpm}
		"""

rule getHe:
	output: outDir + "/he/raw/{accession}_1.fastq.gz", outDir + "/he/raw/{accession}_2.fastq.gz"
	resources: mem=5, time=60*25*3
	shell:
		"""
		module load sra-tools/2.8.2-1
		fastq-dump --gzip --split-files  --outdir {outDir}/he/raw {wildcards.accession}
		"""

rule trimHe:
	input: 
		read1 = outDir + "/he/raw/{accession}_1.fastq.gz",
		read2 = outDir + "/he/raw/{accession}_2.fastq.gz"
	output:
		read1 = outDir + "/he/reads_trimmed/{accession}_R1.fastq.gz",
		read2 = outDir + "/he/reads_trimmed/{accession}_R2.fastq.gz",
	resources: mem=10, time=60*24
	threads: 6
	shell:
		"""
		module load adapterremoval/2.2.2
		AdapterRemoval --trimwindows 5 --minquality 20 --file1 {input.read1} --file2 {input.read2} --output1 {output.read1} --output2 {output.read2} --basename {outDir}/he/reads_trimmed/{wildcards.accession} --gzip --threads {threads}
		"""

rule rsemHe:
	input: rules.trimHe.output
	output: outDir + "/he/rsem/{accession}.genes.results"
	resources: mem=80, time=60*4
	threads: 12
	shell:
		"""
		module load rsem bowtie
		cd {outDir}/he/rsem
		rsem-calculate-expression -p {threads} --paired-end <(zcat {input[0]}) <(zcat {input[1]}) {outDir}/../ref/hbee {wildcards.accession}
		"""

rule collectRsemHe:
	input: expand(outDir + "/he/rsem/{accession}.genes.results", accession = HE.keys())
	output: 
		counts = outDir + "/he/rsem/he_genes_counts.csv", 
		tpm = outDir + "/he/rsem/he_tpm.csv"
	shell:
		"""
		python scripts/collect_counts.py genes {outDir}/he/rsem/ > {output.counts}
		python scripts/collect_tpm.py genes {outDir}/he/rsem/ > {output.tpm}
		"""


rule trim:
	input: outDir + "/reads/{sample}.fastq.gz",
	output: outDir + "/reads_trimmed/{sample}.fastq.gz"
	resources: mem=5, time=60
	threads: 1
	shell:
		"""
		cutadapt -q 10 --trim-n -m 35 -f fastq \
			-a PolyA=AATTGCAGTGGTATCAACGCAGAGCGGCCGC \
			-a TS=AAGCAGTGGTATCAACGCAGAGTACATGGGG \
			-g PolyArc=GCGGCCGCTCTGCGTTGATACCACTGCAATT \
			-g TSrc=CCCCATGTACTCTGCGTTGATACCACTGCTT \
			-g nextera=CTGTCTCTTATACACATCT \
			-o {output} {input}
		"""

rule kallisto:
	input: rules.trim.output
	output: outDir + "/kallisto/{sample}/abundance.tsv"
	resources: mem=5, time=60
	threads: 4
	shell:
		"""
		module load kallisto
		kallisto quant -t {threads} --single -s 100 -l 350 -b 100 -i {kallistoIndex} -o {outDir}/kallisto/{wildcards.sample} {input}
		"""

rule rsem:
	input: rules.trim.output
	output: outDir + "/rsem/{sample}.genes.results"
	resources: mem=10, time=60*2
	threads: 12
	shell:
		"""
		module load rsem bowtie
		cd {outDir}/rsem
		rsem-calculate-expression -p {threads} <(zcat {input}) {outDir}/../ref/hbee {wildcards.sample}
		"""

rule collectRsem:
	input: expand(outDir + "/rsem/{sample}.genes.results", sample = SAMPLES)
	output: outDir + "/rsem/genes.results"
	shell:
		"""
		module load rsem
		rsem-generate-data-matrix {input} > {output}
		"""
rule rsemEbseq:
	input: rules.collectRsem.output
	output:
		ebseq = outDir + "/rsem/genes.ebseq",
		ebseqFDR = outDir + "/rsem/genes.ebseq.fdr"
	shell:
		"""
		module load rsem
		rsem-run-ebseq {input} 10,10 {output.ebseq}
		rsem-control-fdr {output.ebseq} 0.05 {output.ebseqFDR}
		"""
'''Snakefile for RNA pol II ChIP-seq class analysis'''

from datetime import date
today = date.today().isoformat()
workdir: 'results-polII-chip-seq'

#config file supplies sample information and should supply defined
#variables as well but those arenot able to be found by snakemake
#unless defined in Snakefile as shown below

configfile: "/vol3/home/ancarr/projects/pol2_chipseq/config.yaml"


SAMPLES = config['SAMPLES']

FASTA = '/vol3/home/ancarr/data-sets/genome/hg19.fa.gz'
GENES = '/vol3/home/ancarr/data-sets/genome/knownGene.hg19.bed.gz'
CHROMS = '/vol3/home/ancarr/data-sets/genome/hg19.genome'
BOWTIEINDEX = '/vol3/home/jhessel/ref/genomes/hg19/hg19'
WINDOWS = '/vol3/home/ancarr/data-sets/genome/tss_windows.bed'

#final sample type is bed file with signal coverage maps to 2 Kb regions
#around the TSS

rule all:
  input: expand("bedfiles/{samples}/coverage.bed", samples=SAMPLES)


#maps signal from the bedgraph files to 2kb "windows" around the TSS
#use a bedfile I made with the window regions defined for each gene

rule get_coverage:
  input: 
    windows = WINDOWS,
    bedgraph = "bedgraphs/{samples}.bedgraph"
  output:
    coverage = "bedfiles/{samples}/coverage.bed"
  params:
    job_name = 'coverage',
    memory =  "select[mem>4] rusage[mem=4]"
  log:
    'logs/{samples}/coverage.log'
  shell:
    """
    bedtools map -a {input.windows} \
      -b {input.bedgraph} \
      -c 4 -o mean -null 0 > {output.coverage}
    """
 
#generates bedgraph files from bam files to map signal to the genome

rule make_bedgraphs:
  input:
    "bowtie/{samples}.bam"
  output:
    "bedgraphs/{samples}.bedgraph"
  params:
    job_name = 'bedgraph',
    memory =  "select[mem>4] rusage[mem=4]"
  log:
    'logs/{samples}/bedgraph.log'
  shell:
    """
    bedtools genomecov \
      -ibam {input} \
      -bg \
      -g {CHROMS} > {output}
    """

#aligns sample fastq files using bowtie2

rule bowtie_mapping:
  input:
    fq = "/vol3/home/ancarr/data/chip-seq/160115_pol2_chip/{samples}.fastq.gz"
  output:
    "bowtie/{samples}.bam"
  params:
    idx = "/vol3/home/jhessel/ref/genomes/hg19/hg19",
    memory =  "select[mem>40] rusage[mem=40]",
    job_name = 'align'
  threads:
    8
  log:
    "logs/{samples}/align.log"
  shell:
    """
    bowtie2 \
      -t {threads} \
      -x {params.idx} \
      -U {input.fq} \
      -S {output}.tmp
    
    samtools sort {output}.tmp > {output}
    samtools index {output}
 
    """

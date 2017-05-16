'''Snakefile for RNA pol II ChIP-seq class analysis'''

from datetime import date
today = date.today().isoformat()
workdir: 'results-polII-chip-seq'

configfile: "/vol3/home/ancarr/projects/pol2_chipseq/config.yaml"

#FASTA = config['FASTA']
#CHROMS = config['CHROMS']
#GENES = config['GENES']
SAMPLES = config['SAMPLES']

FASTA = '/vol3/home/ancarr/data-sets/genome/hg19.fa.gz'
GENES = '/vol3/home/ancarr/data-sets/genome/knownGene.hg19.bed.gz'
CHROMS = '/vol3/home/ancarr/data-sets/genome/hg19.genome'
BOWTIEINDEX = '/vol3/home/jhessel/ref/genomes/hg19/hg19'
WINDOWS = '/vol3/home/ancarr/data-sets/genome/tss_windows.bed'


rule all:
  input: expand("bedfiles/{samples}/coverage.bed", samples=SAMPLES)

rule get_coverage:
  input: 
    windows = WINDOWS
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

rule prepare_windows:
  input:
    tss = TSS,
    genes = GENES,
    chroms = CHROMS
  output:
    "bedfiles/tss_windows.bed"
  params:
    job_name = 'windows',
    memory =  "select[mem>4] rusage[mem=4]"
  log:
    'logs/windows.log'
  shell:
    "{input.tss} {input.genes} {input.chroms} {output}"
 
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

   # FASTA 
  #output:
   # "indexes/bowtie_idx/hg19.1.bt2"
  #params:
   # job_name = 'indexes',
    #output_name = "indexes/bowtie_idx/hg19",
    #memory =  "select[mem>30] rusage[mem=30] "
  #threads:
   # 8
  #log:
   # "logs/bowtie2_idx.log"
  #shell:
   # """
   # gunzip {input} -c > {input}.tmp
    #bowtie2-build {input}.tmp {params.output_name}
   # """

rule bowtie_mapping:
  input:
    fq = "/vol3/home/ancarr/data/chip-seq/160115_pol2_chip/{samples}.fastq.gz",
#    idx = "/vol3/home/jhessel/ref/genomes/hg19/hg19"
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
      -x {params.idx} \
      -U {input.fq} \
      -S {output}.tmp
    
    samtools sort {output}.tmp > {output}
    samtools index {output}
 
    """

#! /usr/bin/env bash

#BSUB -J polII
#BSUB -o results-polII-chip-seq/logs/polII_log.%J.out
#BSUB -e results-polII-chip-seq/logs/polII_log.%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"

set -o nounset -o pipefail -o errexit -x

workdir="results-polII-chip-seq"
if [[ ! -d $workdir/log ]]; then
     mkdir -p $workdir/log
fi

args=' -C 0 -q normal -n {threads} -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory}"'
snakemake --drmaa "$args" --jobs 8 --configfile config.yaml






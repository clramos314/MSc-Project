import pandas as pd
import os
from snakemake.utils import validate

#singularity: "docker://continuumio/miniconda3:4.6.14"

#report: "report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

include: "rules/common.smk"

##### Target rules #####

#rule all:
#    input:
#        f"{OUTDIR}/filtered/all.vep.vcf"

##rule filter_vep:
##    input:
##        f"{OUTDIR}/filtered/all.vep.vcf"
##    output:
##        f"{OUTDIR}/refiltered/all.vep.refiltered.vcf"
##    conda:"envs/vep.yaml"
##    threads: get_resource("filter_vep","threads")
##    resources:
##        mem = get_resource("filter_vep","mem"),
##        walltime = get_resource("filter_vep","walltime")
##    log:
##        f"{LOGDIR}/filter_vep/filter_vep.log"
##    shell: """
##        ./filter_vep -i {input} -filter "IMPACT is HIGH or IMPACT is MODERATE and VARIANT_CLASS is SNV and (Consequence is transcript_ablation or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is stop_lost or Consequence is start_lost or Consequence is transcript_amplification or Consequence is inframe_insertion or Consequence is inframe_deletion or Consequence is missense_variant or Consequence is protein_altering_variant or Consequence is splice_region_variant or Consequence is incomplete_terminal_codon_variant or Consequence is stop_retained_variant) and SIFT is deleterious and Polyphen is possibly_damaging" --only_matched --force_overwrite -o {output}
##    """

rule pyclone_vi:
    input:
        f"../data/pyclone-vi/input/TCGA-38-4629.IMPACT_is_HIGH.vcf"
    output:
        fit = f"../data/pyclone-vi/output/TCGA-38-4629.IMPACT_is_HIGH.refiltered.h5",
        result = f"../data/pyclone-vi/output/TCGA-38-4629.IMPACT_is_HIGH.tsv"
    conda:"envs/pyclone-vi.yaml"
    threads: get_resource("pyclone","threads")
    resources:
        mem = get_resource("pyclone","mem"),
        walltime = get_resource("pyclone","walltime")
    shell:"""
        pyclone-vi fit -i {input} -o {output.fit} -c 40 -d beta-binomial -r 10 &&
        pyclone-vi write-results-file -i {output.fit} -o {output.result}
    """

##### Modules #####

#include: "rules/refiltering.smk"
#include: "rules/pyclone-vi.smk"

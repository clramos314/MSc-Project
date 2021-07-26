rule filter_vep:
    input:
        f"{OUTDIR}/filtered/all.vep.vcf"
    output:
        f"{OUTDIR}/refiltered/all.vep.refiltered.vcf"
    conda:"envs/vep.yaml"
    threads: get_resource("filter_vep","threads")
    resources:
        mem = get_resource("filter_vep","mem"),
        walltime = get_resource("filter_vep","walltime")
    log:
        f"{LOGDIR}/filter_vep/filter_vep.log"
    shell: """
        echo $PWD && filter_vep -i {input} -filter "IMPACT is HIGH or IMPACT is MODERATE and VARIANT_CLASS is SNV and (Consequence is transcript_ablation or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is stop_lost or Consequence is start_lost or Consequence is transcript_amplification or Consequence is inframe_insertion or Consequence is inframe_deletion or Consequence is missense_variant or Consequence is protein_altering_variant or Consequence is splice_region_variant or Consequence is incomplete_terminal_codon_variant or Consequence is stop_retained_variant) and SIFT is deleterious and Polyphen is possibly_damaging" --only_matched --force_overwrite -o {output}
    """

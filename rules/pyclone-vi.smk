rule pyclone_vi:
    input:
        f"../data/pyclone-vi/input/allSRR.tsv"
    output:
        fit = f"../data/pyclone-vi/output/allSRR.h5",
        result = f"../data/pyclone-vi/input/allSRR_pyclone_vi_output.tsv"
    conda:"envs/pyclone-vi.yaml"
    threads: get_resource("pyclone","threads")
    resources:
        mem = get_resource("pyclone","mem"),
        walltime = get_resource("pyclone","walltime")
    shell:"""
        pyclone-vi fit -i {input} -o {output.fit} -c 40 -d beta-binomial -r 10 &&
        pyclone-vi write-results-file -i {output.fit} -o {output.result}
    """

rule donor_align:
    input:
        "results/by_donor/{donor}.fa"
    output:
        "results/by_donor/{donor}.aln"
    conda:
        "../envs/myenv.yaml"
    shell:
        "mafft {input} > {output}"

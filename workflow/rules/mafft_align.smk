rule mafft_align:
    input:
        "results/annot.fa"
    output:
        "results/sequences.aln"
    conda:
        "../envs/myenv.yaml"
    shell:
        "mafft {input} > {output}"

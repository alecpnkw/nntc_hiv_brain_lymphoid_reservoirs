rule fasttree:
    input:
        "results/sequences.aln"
    output:
        "results/sequences.newick"
    conda:
        "../envs/myenv.yaml"
    shell:
        "fasttree -nt {input} > {output}"

rule functional_filter:
   input:
      "resources/sequences.fa"
   output:
      "results/annot.fa",
      "results/intact.fa",
      "results/defective.fa"
   params:
      max_first_stop = 750,
      min_aa_length = 800,
      frame =  3
   script:
      "../scripts/functional_filter.jl"

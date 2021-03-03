# Python implementation from https://gist.github.com/philiptzou/8d6c7c61d2242f730a1a2f87ba9a2a72
rule hypermut:
   input: "results/hypermut_in_{donor}.aln"
   output: temporary("results/hypermut_{donor}.tsv")
   script: "../scripts/hypermut2.py"

# HYPERMUT requires a consensus or refseq as the first record
rule hypermut_input:
   input: "results/by_donor/{donor}.aln"
   output: temporary("results/hypermut_in_{donor}.aln")
   run:
      from Bio import SeqIO, Seq, SeqRecord
      from collections import Counter
      import numpy as np

      # for simple consensus, ignoring ties
      def mode(collection):
         c = Counter(collection)
         return c.most_common(1)[0][0]

      recs = list(SeqIO.parse(str(input), "fasta"))
      seqs = [list(r.seq.upper()) for r in recs]
      seq_array = np.array(seqs)

      cons = []
      for i in range(seq_array.shape[1]):
         cons += list(mode(seq_array[:,i]))
      cons_record = SeqRecord.SeqRecord(Seq.Seq("".join(cons)), id = "CONSENSUS", description = "")
      SeqIO.write([cons_record] + recs, output[0], "fasta")

# Collect into single table
rule gather_hypermut:
   input: expand("results/hypermut_{donor}.tsv", donor = DONORS)
   output: "results/hypermut.tsv"
   run:
      import pandas as pd
      dfs = [pd.read_table(f, sep = '\t') for f in input]
      pd.concat(dfs).to_csv(output[0], sep = '\t', index = False)

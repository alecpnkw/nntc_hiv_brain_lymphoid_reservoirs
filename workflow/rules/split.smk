rule split:
    input:
        "results/annot.fa"
    output:
        temporary("results/by_donor/{donor}.fa")
    run:
        from Bio import SeqIO
        import sys
        import re
        import os

        out_recs = []
        out_dir = os.path.dirname(str(output))

        with open(str(input), 'r') as io:
            for r in SeqIO.parse(io, "fasta"):
                donor = re.match("NNTC_[0-9]{2}", r.id).group()
                if donor == wildcards.donor:
                    out_recs.append(r)

        SeqIO.write(out_recs, str(output),"fasta")

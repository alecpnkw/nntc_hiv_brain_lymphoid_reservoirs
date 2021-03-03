using NextGenSeqUtils: read_fasta_with_names,
write_fasta, translate_to_aa

function robust_translate(s)
    s = replace(s,"-"=>"")
    s = s[1:3*Int64(floor(length(s)/3))]
    translate_to_aa(s)
end

function find_first_stop(AAseq)
    s = collect(AAseq)
    if '*' in s
        return findfirst(s .== '*')
    else
        return length(AAseq)
    end
end

seqnames, seqs = read_fasta_with_names(snakemake.input[1]);
res = []

for (i,s) in enumerate(seqs)
    early_stop, deletion = false, false
    annot_name = seqnames[i]
    aa_seq = robust_translate.(s[snakemake.params["frame"]:end])
    first_stop = find_first_stop(aa_seq)
    if first_stop < snakemake.params["max_first_stop"]
        early_stop = true
        annot_name = annot_name * " early_stop:AA$(first_stop)"
    end
    if length(aa_seq) < snakemake.params["min_aa_length"]
        deletion = true
        annot_name = annot_name * " deletion:$(length(aa_seq))AAs"
    end
    push!(res, (s, aa_seq, annot_name, early_stop, deletion))
end

#all annotated
write_fasta(snakemake.output[1], getindex.(res, 1); names = getindex.(res, 3));

#intact
is_intact = .!(getindex.(res, 4) .| getindex.(res, 4))
write_fasta(
    snakemake.output[2],
    getindex.(res, 1)[is_intact];
    names = getindex.(res, 3)[is_intact]
)

#defective
write_fasta(
    snakemake.output[3],
    getindex.(res, 1)[.!is_intact];
    names = getindex.(res, 3)[.!is_intact]
)

function infernal_align_fasta_to_cm(fasta_path::AbstractString, cm_path::AbstractString)
    # trimmed (no inserts) aligned fasta
    aln = Infernal.cmalign(cm_path, fasta_path; matchonly=true, outformat="AFA")

    # these are already aligned and without inserts
    sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(aln.out)))
    @assert allequal(map(length, sequences))

    return LongRNA{4}.(sequences)
end

function infernal_score_sequences(cm_path::AbstractString, sequences::AbstractVector; informat=nothing, notrunc=true, glob=true)
    mktemp() do fasta_path, io
        FASTX.FASTA.Writer(io) do writer
            for (i, sequence) = enumerate(sequences)
                record = FASTX.FASTA.Record(string(i), sequence)
                write(writer, record)
            end
        end
        result = Infernal.cmalign(cm_path, fasta_path; glob, informat, notrunc)
        return Infernal.cmalign_parse_sfile(result.sfile)
    end
end

function infernal_cm_emit_sequences(cm_path::AbstractString; N::Int=5000, inserts::Bool=false)
    emitted_sequences_afa = Infernal.cmemit(cm_path; N, aligned=true, outformat="AFA")
    emitted_sequences_with_inserts = FASTX.sequence.(FASTX.FASTA.Reader(open(emitted_sequences_afa.out)))
    if inserts
        return LongRNA{4}.(emitted_sequences_with_inserts)
    else
        return LongRNA{4}.(filter(!=('.'), filter(!islowercase, seq)) for seq = emitted_sequences_with_inserts)
    end
end

function rfam_aligned_path(id::String)
    return joinpath(get_data_directory(), "infernal", "$id.full.noinserts.afa")
end

function rfam_reweights_path(id::String)
    return joinpath(get_data_directory(), "infernal", "$id.full.noinserts.afa.wts")
end

function rfam_seed_stk_path(id::String)
    joinpath(get_data_directory(), "infernal", "SEED", "$id.stk")
end

function rfam_load(id::String)
    path = rfam_aligned_path(id)
    str_sequences = open(FASTX.FASTA.Reader, path) do reader
        LongRNA{4}.(FASTX.sequence.(reader))
    end
    return LongRNA{4}.(str_sequences)
end

"""
    reweights(v)

Computes the reweighting vector for the MSA.
"""
function reweights(
    v::BitArray{3};
    verbose::Bool = false,
    θ::Real = 0.2 # similarity threshold (as fraction of sequence length)
)
    q = size(v, 1) # number of letters, should be 5
    wts, Beff = DCAUtils.compute_weights(potts(v), q, θ; verbose)
    return wts
end

function rfam_save_reweights(id::String)
    v = onehot(rfam_load(id))
    wts = reweights(v)
    path = rfam_reweights_path(id)
    writedlm(path, wts)
    return path
end

function rfam_load_reweights(id::String)
    path = rfam_reweights_path(id)
    return vec(readdlm(path, Float))
end

"""
    stockholm_ss(path)

Extract SS string in WUSS format from a Stockholm file.
"""
function stockholm_ss(path::AbstractString)
    ss = ""
    for line in eachline(path)
        if startswith(line, "#=GC SS_cons")
            words = split(line)
            @assert length(words) == 3
            ss *= last(words)
        end
    end
    return ss
end

"""
    rfam_ss(id; inserts=false)

Fetch secondary structure WUSS string from Rfam Id. Pass `inserts=true` to include insertions.
"""
function rfam_ss(id::AbstractString; inserts::Bool=false)
    stk = Infernal.esl_afetch(Rfam.seed(), id)
    ss = stockholm_ss(stk.out)
    if inserts
        return ss # include insertions
    else
        return replace(ss, '.' => "") # remove insertions
    end
end

"""
    rfam_cmenone(id)

Rfam covariance models are heavily regularized to enable them to capture distant homologues.
This means that their statistics will be noisy and a bit far from the natural sequences.
Here we retrain a new CM model on the full alignment (including hits), and without the regularization.
The resulting model should match single-site statistics of the alignment perfectly.
"""
function rfam_cmenone(id::AbstractString)
    # CM model from Rfam (this has the noisy floor!)
    cm = Infernal.cmfetch(Rfam.cm(), id)
    # aligned hits, used to train a new noiseless CM model (in Stockholm format, without inserts!)
    hits_stk = Infernal.cmalign(cm.out, Rfam.fasta_file(id); matchonly=true)
    # fit new CM model using full alignment (without inserts), and without entropic noise
    return Infernal.cmbuild(hits_stk.out; enone=true)
end

function rfam_RF00162_hits()
    fasta = Rfam.fasta_file("RF00162")
    cm = Infernal.cmfetch(Rfam.cm(), "RF00162").out

    sequences = infernal_align_fasta_to_cm(fasta, cm)
    @assert only(unique(map(length, sequences))) == 108
    return sequences
end

function rfam_RF00162_seed()
    RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")
    RF00162_seed_match_cols = findall(≠('.'), stockholm_ss(RF00162_seed_stk.out))

     # WARNING: this has inserts marked as '-'
    RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM")

    RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))
    RF00162_seed_seqs_noinserts = [FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records];

    @assert only(unique(length.(RF00162_seed_seqs_noinserts))) == 108

    return LongRNA{4}.(RF00162_seed_seqs_noinserts)
end

function rfam_RF00162_rfam_cm()
    return Infernal.cmfetch(Rfam.cm(), "RF00162")
end

function rfam_RF00162_denoised_cm()
    return rfam_cmenone("RF00162")
end

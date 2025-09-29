#=
Define the order in which amino-acids and nucleotides are encoded as integers.
=#

function alphabet(::Type{<:LongRNA}; gap::Bool = true)
    NTs = rna"ACGU"
    if gap
        return NTs * rna"-"
    else
        return NTs
    end
end

function alphabet(s::LongRNA; gap::Bool = true)
    return alphabet(typeof(s); gap)
end

#=
Encodes sequences as Potts arrays, where the integer entry indicates the sequence letter.
=#

function potts(X::BitMatrix)
    return vec(Int8.(first.(Tuple.(argmax(X; dims=1)))))
end

function potts(X::BitArray{3})
    return reshape(Int8.(first.(Tuple.(argmax(X; dims=1)))), size(X,2), size(X,3))
end

function potts(s::Union{LongRNA{4}, AbstractVector{<:LongRNA}})
    potts(onehot(s))
end

function rnaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongRNA{4}([alphabet(LongRNA; gap)[i] for i in P])
end

function rnaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [rnaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function onehot(seq::LongRNA{4}; gap::Bool=true)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(alphabet(seq; gap))
end

function onehot(seqs::AbstractVector{<:LongRNA}; gap::Bool=true)
    L = only(unique(length.(seqs))) # all sequences must have same length
    return reshape(mapreduce(s -> onehot(s; gap), hcat, seqs), :, L, length(seqs))
end

function rnaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 5
        gap = true
    elseif size(X, 1) == 4
        gap = false
    else
        error("Expected 4 or 5 rows; got $(size(X, 1))")
    end

    return rnaseq(potts(X); gap)
end

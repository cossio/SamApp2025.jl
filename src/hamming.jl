"""
    hamming(seq1, seq2)

Hamming distance between two sequences.
"""
function hamming(seq1::BitMatrix, seq2::BitMatrix)
    @assert size(seq1) == size(seq2)
    return sum(any(seq1 .!= seq2, dims=1))
end

"""
    hamming(seqs1, seqs2)

Matrix of Hamming inter-distances from `seqs1` to `seqs2`.
"""
function hamming(s1::BitArray{3}, s2::BitArray{3})
    @assert size(s1, 1) == size(s2, 1)
    @assert size(s1, 2) == size(s2, 2)
    return [hamming(s1[:,:,n], s2[:,:,m]) for n in axes(s1, 3), m in axes(s2, 3)]
end

"""
    hamming(seqs)

Matrix of Hamming intra-distances within `seqs`.
"""
function hamming(seqs::BitArray{3})
    d = zeros(Int, size(seqs, 3), size(seqs, 3))
    for n in axes(seqs, 3), m in 1:(n - 1)
        d[m, n] = d[n, m] = hamming(seqs[:,:,m], seqs[:,:,n])
    end
    return d
end

function hamming(seqs::BitArray{3}, seq::BitMatrix)
    @assert (size(seqs, 1), size(seqs, 2)) == size(seq)
    return [hamming(seqs[:,:,n], seq) for n in axes(seqs, 3)]
end

hamming(seq::BitMatrix, seqs::BitArray{3}) = hamming(seqs, seq)'

function hamming(s1::LongRNA, s2::LongRNA)
    @assert length(s1) == length(s2)
    return StringDistances.Hamming()(string(s1), string(s2))
end

hamming(seqs1::AbstractVector{<:LongSequence}, s2::LongRNA) = [hamming(s1, s2) for s1 in seqs1]
hamming(s1::LongRNA, seqs2::AbstractVector) = hamming(seqs2, s1)'

function hamming(seqs1::AbstractVector{<:LongSequence}, seqs2::AbstractVector{<:LongSequence})
    return [hamming(s1, s2) for s1 in seqs1, s2 in seqs2]
end

function hamming(seqs::AbstractVector{<:LongSequence})
    d = zeros(Int, length(seqs), length(seqs))
    for (n, s1) in enumerate(seqs)
        for (m, s2) in enumerate(seqs)
            if m ≥ n
                break
            else
                d[m, n] = d[n, m] = hamming(s1, s2)
            end
        end
    end
    return d
end

levenshtein(s1::LongRNA, s2::LongRNA) = Levenshtein()(string(s1), string(s2))
levenshtein(seqs1::AbstractVector{<:LongSequence}, s2::LongRNA) = [levenshtein(s1, s2) for s1 in seqs1]
levenshtein(s1::LongRNA, seqs2::AbstractVector{<:LongSequence}) = levenshtein(seqs2, s1)'

function levenshtein(seqs1::AbstractVector{<:LongSequence}, seqs2::AbstractVector{<:LongSequence})
    return [levenshtein(s1, s2) for s1 in seqs1, s2 in seqs2]
end

function levenshtein(seqs::AbstractVector{<:LongSequence})
    d = zeros(Int, length(seqs), length(seqs))
    for (n, s1) in enumerate(seqs)
        for (m, s2) in enumerate(seqs)
            if m ≥ n
                break
            else
                d[m, n] = d[n, m] = levenshtein(s1, s2)
            end
        end
    end
    return d
end

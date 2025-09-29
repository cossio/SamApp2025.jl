function effective_couplings(rbm::RBM, data::AbstractArray; wts=nothing)
    @assert size(data) == (size(rbm.visible)..., size(data)[end])
    ν = batchmean(rbm.hidden, var_h_from_v(rbm, data); wts)
    w = flat_w(rbm)
    J = w * Diagonal(vec(ν)) * w'
    return reshape(J, size(rbm.visible)..., size(rbm.visible)...)
end

"""
    effective_contacts(rbm, data)

Compute the DCA-like effective contacts, for an RBM and a MSA dataset. See
Tubiana et al Elife'2019 .
"""
function effective_contacts(rbm::RBM{<:Potts}, data::AbstractArray; wts=nothing)
    J = effective_couplings(rbm, data; wts)
    A = sum(Base.Fix2(^, 2), J; dims=(1, ndims(rbm.visible) + 1))
    L = prod(size(rbm.visible)[2:end])
    apc = average_product_correction(reshape(A, L, L))
    return reshape(apc, size(rbm.visible)[2:end]..., size(rbm.visible)[2:end]...)
end

"""
    average_product_correction(matrix)

Compute the average product correction (APC) for a matrix. See Dunn et al 2008.
"""
function average_product_correction(A::AbstractMatrix)
    @assert size(A, 1) == size(A, 2)
    return A .- sum(A; dims=1) .* sum(A; dims=2) ./ sum(A)
end

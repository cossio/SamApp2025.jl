"""
    RF00162_sites_paired()

Returns the lists of paired, unpaired, and pseudoknot sites in the consensus secondary
structure of RF00162.
"""
function RF00162_sites_paired()
    wuss = RF00162_wuss(; insertions=false)
    bps = findall(∈("()<>"), wuss)
    nps = findall(∉("()<>Aa"), wuss)
    pks = findall(∈("Aa"), wuss)
    return (; bps, nps, pks)
end

"""
    RF00162_wuss(; insertions=false)

Fetch the WUSS secondary structure of RF00162 from Rfam.
"""
function RF00162_wuss(; insertions=false)
    wuss_full = stockholm_ss(Infernal.esl_afetch(Rfam.seed(), "RF00162").out)
    if insertions
        return wuss_full
    else
        # remove insertions
        return filter(!=('.'), wuss_full)
    end
end

"""
    RF00162_sites_annotated_secondary_structure()

Returns the lists of sites in the consensus secondary structure of RF00162, in different
structural elements (p1, p2, ..., base pairs, unpaired, pseudoknots, ...).
"""
function RF00162_sites_annotated_secondary_structure()
    wuss = RF00162_wuss(; insertions=false)
    bps = findall(∈("()<>"), wuss)
    nps = findall(∉("()<>Aa"), wuss)
    pk = findall(∈("Aa"), wuss)
    #= Note that the pseudoknot sites are not included in bps nor nps. =#
    p1 = union(1:8, 101:108)
    p2 = union(13:17, 38:42, 21:23, 29:31)
    p3 = union(43:46, 48:53, 61:64, 67:72)
    p4 = union(81:86, 92:97)
    return (; bps, nps, pk, p1, p2, p3, p4)
end

"""
    clean_wuss(wuss; keep_pseudoknots = false)

Replaces all brackets by parentheses, and all unpaired sites by dots.
"""
function clean_wuss(wuss::AbstractString; keep_pseudoknots::Bool = false)
    _wuss = replace(
        wuss,
        r"\(|\[|\{|\<" => '(',
        r"\)|\]|\}|\>" => ')',
        r"\-|\_|\," => '.'
    )
    if keep_pseudoknots
        return _wuss
    else
        return replace(_wuss, r"[A-Z]|[a-z]" => '.')
    end
end

# RNAeval for the pseudoknot
function vienna_pk_binding_energy_rnaeval(seq::String)
    # store here so that we don't have to fetch this every time
    wuss = "((((((((,,,,<<<<<---<<<_AAAA>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<aaaa<<<<<<_____>>>>>>-->))))))))"

    # only Pk is base-paired here
    ss_pk_only = replace(wuss, r"\(|\)|\[|\]|\{|\}|\<|\>|\-|\_|\," => '.', 'A' => '(', 'a' => ')')

    _rnaeval_input, io = mktemp()
    write(io, seq * '\n')
    write(io, ss_pk_only * '\n')
    close(io)

    _rnaeval_stdout = tempname()
    run(pipeline(`$(ViennaRNA_jll.RNAeval()) -vi $_rnaeval_input`; stdout=_rnaeval_stdout))

    _lines = readlines(_rnaeval_stdout)[1:5]
    _vals = [parse(Int, split(l, [':', ' ']; keepempty=false)[end]) for l in _lines]

    return sum(_vals[2:4]) / 100
end

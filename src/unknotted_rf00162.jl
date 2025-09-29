#=
By permuting the columns of the RF00162, we can "unknot" the pseudoknot, so that it
can be modeled in the Covariance Model as part of the secondary structure.
=#

function rfam_RF00162_unknotted_permutation()
    wuss = rfam_ss("RF00162")
    return [
        1:findlast(==('A'), wuss); # first segment up to first branch of pseudoknot
        findfirst(==('a'), wuss):findlast(==('a'), wuss); # second branch of pseudoknot
        (findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1); # what's in between the pseudoknot
        (findlast(==('a'), wuss) + 1):108 # what's after the pseudoknot
    ]
end

function rfam_RF00162_unknotted_hits_stk()
    Rfam_cm = SamApp2025.rfam_RF00162_rfam_cm()
    RF00162_hits_stk = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true);

    wuss = rfam_ss("RF00162")
    perm = rfam_RF00162_unknotted_permutation()

    # build training alignment for the Untangled CM
    RF00162_hits_stk_permuted = tempname()

    # RF00162_hits_stk has the hits in Stockholm format, with one line per sequence
    open(RF00162_hits_stk_permuted, "w") do file
        for (line_index, line) in enumerate(eachline(RF00162_hits_stk.out))
            if startswith(line, "# STOCKHOLM") || startswith(line, "#=GF") || startswith(line, "#=GS") || isempty(line) || line == "//"
                write(file, line, '\n') # these comment lines we just copy
            elseif startswith(line, "#=GR")
                continue # skip GR annotations
            elseif startswith(line, "#=GC SS_cons") # consensus secondary structure
                _ss = line[41:end]
                @assert _ss  == "((((((((,,,,<<<<<---<<<_____>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<----<<<<<<_____>>>>>>-->))))))))"
                @assert wuss == "((((((((,,,,<<<<<---<<<_AAAA>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<aaaa<<<<<<_____>>>>>>-->))))))))"
                @assert _ss == replace(wuss, 'A' => '_', 'a' => '-')
                # add parenthesis for pseudoknot
                _ss = _ss[1:(findfirst(==('A'), wuss) - 1)] * "((((" * _ss[(findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1)] * "))))" * _ss[(findlast(==('a'), wuss) + 1):end]
                @assert length(_ss) == 108
                write(file, line[1:40] * _ss[perm], '\n') # write untangled (permuted) secondary structure
            elseif startswith(line, "#=GC RF") # consensus sequence
                @assert line[41:end] == "cucUuAUcaAGAGgGGcgGAGGGAcuGGCCCuaUGAAgCCcCgGCAACCccccauaauaaggggAaGGUGCcAAuuCCugCcggccauuaaggccgGaaagAUaAgag"
                write(file, line[1:40] * line[41:end][perm], '\n')
            else # this is a sequence line
                @assert length(line[41:end]) == 108
                write(file, line[1:40] * line[41:end][perm], '\n')
            end
        end
    end

    return RF00162_hits_stk_permuted
end

function rfam_RF00162_unknotted_cm()
    RF00162_hits_stk_permuted = rfam_RF00162_unknotted_hits_stk()
    return Infernal.cmbuild(RF00162_hits_stk_permuted; enone=true)
end

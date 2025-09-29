import CSV
import FASTX
import Infernal
import Rfam
import SamApp2025
using DataFrames: DataFrame
using Test: @test

natural_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")
natural_seed_afa = Infernal.esl_reformat("afa", natural_seed_stk.out)
natural_seed_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_seed_afa.out)))
seed_ss = SamApp2025.stockholm_ss(natural_seed_stk.out)

natural_seed_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_seed_afa.out)))
natural_seed_sequences_afa = String[]
for s in natural_seed_sequences
    push!(natural_seed_sequences_afa, join([f == '.' ? (c == '-' ? '.' : lowercase(c)) : c for (c, f) in zip(s, seed_ss)]))
end

# hits
natural_hits_afa = Infernal.cmalign(
    Infernal.cmfetch(Rfam.cm(), "RF00162").out,
    Rfam.fasta_file("RF00162");
    outformat="afa"
)
natural_hits_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_hits_afa.out)))
natural_hits_sequences_afa = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_hits_afa.out)))

@test allequal(length.(natural_hits_sequences_afa))
@test allequal(length.(natural_seed_sequences_afa))

natural_hits_sequences_matchonly = [filter(!islowercase, replace(seq, '.' => "")) for seq in natural_hits_sequences_afa]
natural_seed_sequences_matchonly = [filter(!islowercase, replace(seq, '.' => "")) for seq in natural_seed_sequences_afa]

@test only(unique(length.(natural_hits_sequences_matchonly))) == 108
@test only(unique(length.(natural_seed_sequences_matchonly))) == 108

# this is the thing we want to test!
pos_align = SamApp2025.shape_positions_alignment_2022()

aptamers_df = SamApp2025.probed_aptamers_table_20221027()

@test sum(aptamers_df.source .== "RF00162_seed70") == 151
@test sum(aptamers_df.source .== "RF00162_full30") == 55

# Note that sequences in `aptamers_df` are DNA, so we want to replace `T` with `U`
aptamers_df.sequence_rna = replace.(aptamers_df.sequence, 'T' => 'U')

# The idea is that pos_align.natural[n,j] will give the position of match column 'j'
# in the probed aptamer 'n' sequence (or `missing`, if this column is deleted (gap) for this aptamer!
# pos_align.synthetic is the same for synthetic (which are easier because there are no inserts)
@test size(pos_align.natural, 2) == size(pos_align.synthetic, 2) == 108
@test size(pos_align.natural, 1) == sum(aptamers_df.source .== "RF00162_seed70") + sum(aptamers_df.source .== "RF00162_full30")
@test size(pos_align.synthetic, 1) == 100

full_aptamers_index_in_rfam = indexin(aptamers_df.id[aptamers_df.source .== "RF00162_full30"], natural_hits_ids)
seed_aptamers_index_in_rfam = indexin(aptamers_df.id[aptamers_df.source .== "RF00162_seed70"], natural_seed_ids)

@test count(isnothing, full_aptamers_index_in_rfam) == 5
@test count(isnothing, seed_aptamers_index_in_rfam) == 0
@test aptamers_df.id[aptamers_df.source .== "RF00162_seed70"] ⊆ natural_seed_ids

for n in 1:55 # hits
    @test aptamers_df.source[n] == "RF00162_full30"
    idx = full_aptamers_index_in_rfam[n]
    isnothing(idx) && continue
    for j in 1:108
        if ismissing(pos_align.natural[n,j])
            @test natural_hits_sequences_matchonly[idx][j] == '-'
        else
            @test aptamers_df.sequence_rna[n][ pos_align.natural[n,j] ] == natural_hits_sequences_matchonly[idx][j]
        end
    end
end

for n in 56:206 # seed
    @test aptamers_df.source[n] == "RF00162_seed70"
    idx = seed_aptamers_index_in_rfam[n - 55]
    @test !isnothing(idx)
    for j in 1:108
        if ismissing(pos_align.natural[n,j])
            @test natural_seed_sequences_matchonly[idx][j] == '-'
        else
            @test aptamers_df.sequence_rna[n][ pos_align.natural[n,j] ] == natural_seed_sequences_matchonly[idx][j]
        end
    end
end

@test sum(aptamers_df.source .== "RF00162_seed70") + sum(aptamers_df.source .== "RF00162_full30") == 206
@test sum(aptamers_df.source .== "RF00162_syn") == 100

synthetic_df = CSV.read(SamApp2025.pierre20221107post(:synthetic), DataFrame)
@test only(unique(length.(synthetic_df.sequence))) == 108

for n in 1:100 # synthetic
    seq_not_aligned = aptamers_df.sequence_rna[aptamers_df.source .== "RF00162_syn"][n]
    seq_yes_aligned = synthetic_df.sequence[n]

    @test filter(!=('-'), seq_yes_aligned) == seq_not_aligned

    for j in 1:108
        if ismissing(pos_align.synthetic[n,j])
            @test seq_yes_aligned[j] == '-'
        else
            @test seq_not_aligned[ pos_align.synthetic[n,j] ] == seq_yes_aligned[j]
        end
    end
end

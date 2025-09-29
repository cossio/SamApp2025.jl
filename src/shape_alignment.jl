#=
This file is concerned with finding the alignment (position mapping) between the sequences
in the shape reactivity files (which are gapless) and the aligned sequences in Rfam (which
are the positions I use in the model). This is a first step in order to process the
reactivity profiles.

We here take care of aligning natural and synthetic sequences.
=#

#=
MORE COMMENTS

Our first goal is to map positions in the probed sequences to the alignemd (match) positions of our model.
These sequences have been already aligned, but the SHAPE data files have lost the information about this alignment.
My goal here is to produce, for each probed aptamer, a list of positions that extract the match sites from the sequences in the SHAPE data files.

The SHAPE data files (see below) have the same sequences as the column `aptamers_df.sequence_tag_primer`.
These sequences include a tag+primer postfix (appended in the end), that we do not want to keep.
Having removed this postfix, the remaining prefix part should match the Rfam (or our designed) sequences.
=#

# The idea is that natural_positions_mapping[n,j] will give the position of match column 'j'
# in the probed aptamer 'n' sequence (or `missing`, if this column is deleted (gap) for this aptamer!
# natural_positions_mapping is for the natural sequences.
# Synthetic sequences are easier (there are no inserts), and I will do a similar array below.

function shape_positions_alignment_2022()
    # load Rfam family
    natural_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")
    natural_seed_afa = Infernal.esl_reformat("afa", natural_seed_stk.out)
    natural_seed_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_seed_afa.out)))
    seed_ss = stockholm_ss(natural_seed_stk.out)

    # these guys don't distinguish insertions, have all gaps as just '-'
    natural_seed_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_seed_afa.out)))

    # therefore .... indicate insertions with '.' and lowercase in natural seed sequences
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
    natural_hits_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_hits_afa.out)))

    # The following file contains sequences (also + tag and primer) probed in 2022 experiments, and their IDs
    # The column sequence_tag_primer has the same sequence that appears in the SHAPE files below
    aptamers_df = probed_aptamers_table_20221027_rna()

    # Seed IDs match exactly what I get from Rfam
    @assert aptamers_df.id[aptamers_df.source .== "RF00162_seed70"] ⊆ natural_seed_ids

    # For the synthetic sequences ... Aptamer names are sorted, and go from 1 to 100 ... just the ones we sent them
    @assert parse.(Int, aptamers_df.id[aptamers_df.source .== "RF00162_syn"]) == 1:100

    # For the natural sequences, aptamer names (indices) are also sorted simply, with the first being the hits, then the seed sequences
    @assert parse.(Int, [s[7:end] for s = aptamers_df.name[aptamers_df.source .== "RF00162_full30"]]) == 1:55
    @assert parse.(Int, [s[7:end] for s = aptamers_df.name[aptamers_df.source .== "RF00162_seed70"]]) == 56:206
    @assert parse.(Int, [s[7:end] for s = aptamers_df.name[(aptamers_df.source .== "RF00162_full30") .| (aptamers_df.source .== "RF00162_seed70")]]) == 1:206
    @assert sum(aptamers_df.source .== "RF00162_full30") + sum(aptamers_df.source .== "RF00162_seed70") == 206

    # considering only those that match Rfam ids
    N_probed_aptamers_natural = (
        count(∈(natural_seed_ids), aptamers_df.id[aptamers_df.source .== "RF00162_seed70"]) +
        count(∈(natural_hits_ids), aptamers_df.id[aptamers_df.source .== "RF00162_full30"])
    )
    @assert N_probed_aptamers_natural == 201
    # there are 5
    # For hit sequences, there only 5 missmatches .... maybe they used a different Rfam version?
    # Not sure ... They are few so I will just ignore these sequences

    # the sequence_tag_primer is what is reported in the SHAPE files
    # This is just the concatenation of the Aptamer sequence + a tag + a primer
    @assert aptamers_df.sequence_tag_primer == aptamers_df.sequence .* aptamers_df.tag .* aptamers_df.primer

    # index of probed aptamer, in my list of natural seed sequences from Rfam
    aptamer_seed_index = indexin(aptamers_df.id[aptamers_df.source .== "RF00162_seed70"], natural_seed_ids)
    aptamer_seed_index = map(identity, aptamer_seed_index)
    @assert all(!isnothing, aptamer_seed_index) # no nothing here!, all seed match
    aptamer_hits_index = indexin(aptamers_df.id[aptamers_df.source .== "RF00162_full30"], natural_hits_ids)

    # The probed natural seed sequences coincide alright with what I get from Rfam!
    @assert aptamers_df.sequence_rna[aptamers_df.source .== "RF00162_seed70"] == replace.(natural_seed_sequences[aptamer_seed_index], '-' => "")

    # The probed natural hits sequences also coincide alright with what I get from Rfam!
    @assert (
        aptamers_df.sequence_rna[aptamers_df.source .== "RF00162_full30"][map(!isnothing, aptamer_hits_index)] ==
        uppercase.(replace.(natural_hits_sequences[filter(!isnothing, aptamer_hits_index)], '-' => "", '.' => ""))
    )

    # A file I prepared containing the synthetic aptamers we designed, together with the name we use for them in the other files
    # (As well as other info). The sequences in this file have gaps, that we can use to align the sequences given in the other files
    # which are gapless.
    synthetic_df = CSV.read(pierre20221107post(:synthetic), DataFrame)

    # The idea is that natural_positions_mapping[n,j] will give the position of match column 'j'
    # in the probed aptamer 'n' sequence (or `missing`, if this column is deleted (gap) for this aptamer!
    # natural_positions_mapping is for the natural sequences.
    # Synthetic sequences are easier (there are no inserts), and I will do a similar array below.

    # First for hits sequences
    natural_hits_positions_mapping = Array{Union{Missing,Int}}(undef, sum(aptamers_df.source .== "RF00162_full30"), 108);
    natural_hits_positions_mapping .= missing

    for (n, idx) = enumerate(aptamer_hits_index)
        isnothing(idx) && continue # sequence not found in Rfam
        aligned_sequence = replace(natural_hits_sequences[idx], '.' => "")
        seq_pos = match_pos = 0
        for (i,c) = enumerate(aligned_sequence)
            if c == '-'
                match_pos += 1 # deletion
                natural_hits_positions_mapping[n, match_pos] = missing # deletion
            elseif isuppercase(c)
                seq_pos += 1
                match_pos += 1 # match
                natural_hits_positions_mapping[n, match_pos] = seq_pos # deletion
            else
                @assert islowercase(c) # insertion
                seq_pos += 1
            end
        end
    end

    # now same deal ... but for seed sequencse
    natural_seed_positions_mapping = Array{Union{Missing,Int}}(undef, sum(aptamers_df.source .== "RF00162_seed70"), 108)
    natural_seed_positions_mapping .= missing

    for (n, idx) = enumerate(aptamer_seed_index)
        @assert !isnothing(idx) # all seed sequences are found in Rfam!
        aligned_sequence = replace(natural_seed_sequences_afa[idx], '.' => "")
        seq_pos = match_pos = 0
        for (i,c) = enumerate(aligned_sequence)
            if c == '-'
                match_pos += 1 # deletion
                natural_seed_positions_mapping[n, match_pos] = missing
            elseif isuppercase(c)
                seq_pos += 1
                match_pos += 1 # match
                natural_seed_positions_mapping[n, match_pos] = seq_pos
            else
                @assert islowercase(c) # insertion
                seq_pos += 1
            end
        end
    end

    # concat natural sequences position mapping (hits + seed)
    natural_positions_mapping = vcat(natural_hits_positions_mapping, natural_seed_positions_mapping);

    # Same as natural_positions_mapping, but for the 100 synthetic aptamers
    # This is a bit simpler than the natural case, because there are no insertions anywhere
    synthetic_positions_mapping = Array{Union{Missing,Int}}(undef, 100, 108)
    synthetic_positions_mapping .= missing

    for n = 1:100
        aligned_sequence = synthetic_df.sequence[n]
        seq_pos = match_pos = 0
        for (i,c) = enumerate(aligned_sequence)
            if c == '-'
                match_pos += 1 # deletion
                synthetic_positions_mapping[n, match_pos] = missing
            else
                @assert isuppercase(c)
                seq_pos += 1
                match_pos += 1 # match
                synthetic_positions_mapping[n, match_pos] = seq_pos
            end
        end
    end

    pdb_positions_mapping = shape_pdb_positions_mapping_20240731()

    return (; natural = natural_positions_mapping, synthetic = synthetic_positions_mapping, pdb = pdb_positions_mapping)
end


function shape_pdb_positions_mapping_20240731()
	aptamers_df = probed_aptamers_table_20221027_rna()

	natural_hits_afa_1410 = Infernal.cmalign(
	   Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out,
	   Rfam.fasta_file("RF00162"; rfam_version="14.10");
	   outformat="afa"
	)

	natural_hits_dsc_1410 = FASTX.description.(FASTX.FASTA.Reader(open(natural_hits_afa_1410.out)))
	natural_hits_sequences_1410 = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_hits_afa_1410.out)))

	indexes_in_aptamers_df = [
		[i for (i, id) = enumerate(aptamers_df.id) if occursin("2GIS", id)];
		[i for (i, id) = enumerate(aptamers_df.id) if occursin("4KQY", id)]
	]

	indexes_in_hits410 = [
		[i for (i, dsc) = enumerate(natural_hits_dsc_1410) if occursin("2GIS", dsc)];
		[i for (i, dsc) = enumerate(natural_hits_dsc_1410) if occursin("4KQY", dsc)]
	]

	@assert length(indexes_in_aptamers_df) == length(indexes_in_hits410) == 2

	pdb_positions_mapping = Array{Union{Missing,Int}}(undef, 2, 108);
	pdb_positions_mapping .= missing

	for (n, index_in_aptamers_df, index_in_hits410) = zip(1:2, indexes_in_aptamers_df, indexes_in_hits410)
        aligned_sequence = replace(natural_hits_sequences_1410[index_in_hits410], '.' => "")
        seq_pos = match_pos = 0
        for (i, c) = enumerate(aligned_sequence)
            if c == '-'
                match_pos += 1 # deletion
                pdb_positions_mapping[n, match_pos] = missing # deletion
            elseif isuppercase(c)
                seq_pos += 1
                match_pos += 1 # match
                pdb_positions_mapping[n, match_pos] = seq_pos # deletion
            else
                @assert islowercase(c) # insertion
                seq_pos += 1
            end
        end
    end

	return pdb_positions_mapping
end

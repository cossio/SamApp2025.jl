# Data demuxed by Pierre for replicates 4+5 and then for 1,2,3.
# Here we download it as an Artifact from Github gists.
function load_shapemapper_data_pierre_demux_20230929_repl123()
    #shape_dir = joinpath(PIERRE_DEMUX, "SAMAP-rep1+2+3-demux")
    shape_dir = artifact"SAMAP-demux-repl123-2023-09" # only demuxed data for repl 1, 2, 3
    @assert isdir(shape_dir)

    #= Note that repl 1, 2, 3 has only natural data =#

    conditions = filter(startswith("SAMAP"), readdir(shape_dir))
    @assert issorted(conditions)

    aptamers_df = probed_aptamers_table_20221027()
    positions_mapping = shape_positions_alignment_2022()

    # Store reactivities in a 3-dimensional tensor,
    # shape_reactivities[i, n, c] = reactivity of site 'i', of sequence 'n', in condition 'c'
    # `missing` denotes sites without measurement (or deletions in the sequence)
    shape_reactivities_natur = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_natur_err = fill(NaN, 108, 206, length(conditions)) # for natural sequences

    # Reactivity profile after normalization, see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#reactivity-profile-calculation-and-normalization
    shape_reactivities_indivnorm_natur = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_natur_err = fill(NaN, 108, 206, length(conditions))

    shape_raw_reactivities_natur = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_natur_err = fill(NaN, 108, 206, length(conditions))

    shape_M_natur = fill(NaN, 108, 206, length(conditions))
    shape_U_natur = fill(NaN, 108, 206, length(conditions))
    shape_D_natur = fill(NaN, 108, 206, length(conditions))

    # these will contain effective read depths
    # see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#effective-read-depth
    shape_M_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_U_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_D_natur_depth = fill(NaN, 108, 206, length(conditions))

    # aptamer names, consistent with the order in the shape_reactivities tensor
    aptamer_names = ["APSAMN$n" for n in 1:206] # natural only

    # origin of aptamers (RF00162 natural, hits or seed, and the synthetic ones, RBM or Infernal)
    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
    ]
    @assert length(aptamer_names) == length(aptamer_origin)

    missing_files = String[]

    for (c, cond) in enumerate(conditions)
        for n in 1:206 # natural sequences (from RF00162)
            # _profile.txt files (no normalization!)
            prof_file = joinpath(shape_dir, cond, "$(cond)_APSAMN$(n)_profile.txt")
            if isfile(prof_file)
                #@info "Loading $prof_file"
                prof_df = CSV.read(prof_file, DataFrame)
                @assert uppercase(join(prof_df.Sequence)) == replace(aptamers_df.sequence[n], 'T' => 'U')
                for i in 1:108
                    if ismissing(positions_mapping.natural[n,i])
                        continue
                    else
                        j = positions_mapping.natural[n,i]
                    end

                    shape_reactivities_natur[i, n, c] = prof_df.HQ_profile[j]
                    shape_reactivities_natur_err[i, n, c] = prof_df.HQ_stderr[j]

                    shape_raw_reactivities_natur[i, n, c] = prof_df.Reactivity_profile[j]
                    shape_raw_reactivities_natur_err[i, n, c] = prof_df.Std_err[j]

                    shape_M_natur[i,n,c] = prof_df.Modified_rate[j]
                    shape_U_natur[i,n,c] = prof_df.Untreated_rate[j]
                    shape_D_natur[i,n,c] = prof_df.Denatured_rate[j]

                    shape_M_natur_depth[i,n,c] = prof_df.Modified_effective_depth[j]
                    shape_U_natur_depth[i,n,c] = prof_df.Untreated_effective_depth[j]
                    shape_D_natur_depth[i,n,c] = prof_df.Denatured_effective_depth[j]

                    if "Norm_profile" ∈ names(prof_df) # some profiles do not have normalized reactivities because of few reads
                        shape_reactivities_indivnorm_natur[i,n,c] = prof_df.Norm_profile[j]
                        shape_reactivities_indivnorm_natur_err[i,n,c] = prof_df.Norm_stderr[j]
                    end
                end
            else
                @warn "$prof_file not found!"
                push!(missing_files, prof_file)
            end
        end
    end

    # concat things
    # attention! I put natural first, then synthetic, as in aptamers_df
    # NaN denotes missing values (due to alignment, or experimental noise)
    # Repl 1, 2, 3 has no synthetic sequences data!
    shape_reactivities = shape_reactivities_natur
    shape_reactivities_err = shape_reactivities_natur_err

    shape_raw_reactivities = shape_raw_reactivities_natur
    shape_raw_reactivities_err = shape_raw_reactivities_natur_err

    shape_reactivities_indivnorm = shape_reactivities_indivnorm_natur
    shape_reactivities_indivnorm_err = shape_reactivities_indivnorm_natur_err

    # mutation rates
    shape_M = shape_M_natur
    shape_U = shape_U_natur
    shape_D = shape_D_natur

    # read depth
    shape_M_depth = shape_M_natur_depth
    shape_U_depth = shape_U_natur_depth
    shape_D_depth = shape_D_natur_depth

    # standard error of the mutation rates (assuminng Poisson statistics)
    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    # sequence ids
    aptamer_ids = aptamers_df.id[indexin(aptamer_names, aptamers_df.name)]

    # fetch Rfam seed sequences
    RF00162_seed_afa = Infernal.esl_reformat(
        "AFA", Infernal.esl_afetch(Rfam.seed(), "RF00162").out; informat="STOCKHOLM"
    )
    RF00162_seed_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))
    RF00162_seed_seqs_full = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out))))
    RF00162_seed_match_cols = findall(≠('.'), RF00162_wuss(; insertions=true))
    RF00162_seed_seqs = [LongRNA{4}(seq[RF00162_seed_match_cols]) for seq in RF00162_seed_seqs_full]

    # fetch Rfam hits sequences -- already trimmed (no inserts) aligned fasta
    RF00162_hits_afa = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(), "RF00162").out,
        Rfam.fasta_file("RF00162");
        matchonly=true, outformat="AFA"
    )
    RF00162_hits_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out)))
    RF00162_hits_seqs = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))))

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
    @assert length(RF00162_seed_match_cols) == 108

    ## Aligned sequences in indexed in the same ordre as the SHAPE data ... we will fill this array below.
    aligned_sequences = Vector{Union{Missing,String}}(undef, length(aptamer_names))

    # seed sequences (no missmatches, so no need for isnothing check)
    aligned_sequences[aptamer_origin .== "RF00162_seed70"] .= string.(
        RF00162_seed_seqs[indexin(aptamer_ids[aptamer_origin .== "RF00162_seed70"], RF00162_seed_ids)]
    )

    # hits (there are 5 missmatches, therefore we need the `isnothing` check)
    aligned_sequences[aptamer_origin .== "RF00162_full30"] .= [
        isnothing(i) ? missing : string(RF00162_hits_seqs[i]) for i in indexin(aptamer_ids[aptamer_origin .== "RF00162_full30"], RF00162_hits_ids)
    ]

    @assert all(ismissing.(aligned_sequences) .|| (length.(aligned_sequences) .== 108))
    @assert allunique(aptamer_names)
    @assert allunique(aptamer_ids)

    return (; shape_reactivities, shape_reactivities_err,
        shape_raw_reactivities, shape_raw_reactivities_err,
        shape_reactivities_indivnorm, shape_reactivities_indivnorm_err,
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        conditions, aptamer_names, aptamer_origin, aptamer_ids, aligned_sequences,
        missing_files
    )
end

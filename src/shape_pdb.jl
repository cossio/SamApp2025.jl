function load_shapemapper_data_pierre_demux_20240730_with_pdb(; demux::Bool=true)
    shape_dir = shapemapper_data_pierre_demux_20230920_dir(; demux)

    conditions = filter(startswith("SAMAP"), readdir(shape_dir))
    @assert issorted(conditions)

    aptamers_df = probed_aptamers_table_20221027()
    synthetic_df = CSV.read(pierre20221107post(:synthetic), DataFrame)
    positions_mapping = shape_positions_alignment_2022()

    # Store reactivities in a 3-dimensional tensor,
    # shape_reactivities[i, n, c] = reactivity of site 'i', of sequence 'n', in condition 'c'
    # `missing` denotes sites without measurement (or deletions in the sequence)
    shape_reactivities_natur = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences
    shape_reactivities_pdb = fill(NaN, 108, 2, length(conditions)) # for PDB sequences
    shape_reactivities_natur_err = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth_err = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences
    shape_reactivities_pdb_err = fill(NaN, 108, 2, length(conditions)) # for PDB sequences

    # Reactivity profile after normalization, see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#reactivity-profile-calculation-and-normalization
    shape_reactivities_indivnorm_natur = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_pdb = fill(NaN, 108, 2, length(conditions))
    shape_reactivities_indivnorm_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth_err = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_pdb_err = fill(NaN, 108, 2, length(conditions))

    # Raw reactivities (without ShapeMapper quality checks)
    shape_raw_reactivities_natur = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_pdb = fill(NaN, 108, 2, length(conditions))
    shape_raw_reactivities_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth_err = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_pdb_err = fill(NaN, 108, 2, length(conditions))

    shape_M_natur = fill(NaN, 108, 206, length(conditions))
    shape_U_natur = fill(NaN, 108, 206, length(conditions))
    shape_D_natur = fill(NaN, 108, 206, length(conditions))
    shape_M_synth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth = fill(NaN, 108, 100, length(conditions))
    shape_M_pdb = fill(NaN, 108, 2, length(conditions))
    shape_U_pdb = fill(NaN, 108, 2, length(conditions))
    shape_D_pdb = fill(NaN, 108, 2, length(conditions))

    # these will contain effective read depths
    # see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#effective-read-depth
    shape_M_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_U_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_D_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_M_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_M_pdb_depth = fill(NaN, 108, 2, length(conditions))
    shape_U_pdb_depth = fill(NaN, 108, 2, length(conditions))
    shape_D_pdb_depth = fill(NaN, 108, 2, length(conditions))

    # aptamer names, consistent with the order in the shape_reactivities tensor
    aptamer_names = [
        ["APSAMN$n" for n in 1:206]; # natural
        ["APSAMS$n" for n in 1:100]; # synthetic
        ["SAMAP-PDB0", "SAMAP-PDB10"] # PDB
    ]
    # origin of aptamers (RF00162 natural, hits or seed, and the synthetic ones, RBM or Infernal)
    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
        ["RF00162_syn_" * synthetic_df.origin[n] for n in 1:100]; # synthetic (rbm or infernal)
        ["PDB", "PDB"] # from PDB
    ]
    @assert length(aptamer_names) == length(aptamer_origin) == 308

    missing_files = String[]

    for (c, cond) = enumerate(conditions)
        for n = 1:206 # natural sequences (from RF00162)
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

                    if "Norm_profile" ∈ names(prof_df)
                        shape_reactivities_indivnorm_natur[i,n,c] = prof_df.Norm_profile[j]
                        shape_reactivities_indivnorm_natur_err[i,n,c] = prof_df.Norm_stderr[j]
                    end
                end
            else
                @warn "$prof_file not found!"
                push!(missing_files, prof_file)
            end
        end

        for n = 1:100 # synthetic sequences
            # _profile.txt files (no normalization!)
            prof_file = joinpath(shape_dir, cond, "$(cond)_APSAMS$(n)_profile.txt")
            @assert isfile(prof_file) "File $prof_file not found!"
            prof_df = CSV.read(prof_file, DataFrame)
            @assert uppercase(join(prof_df.Sequence)) == replace(aptamers_df.sequence[aptamers_df.source .== "RF00162_syn"][n], 'T' => 'U')
            for i in 1:108
                if ismissing(positions_mapping.synthetic[n,i])
                    continue
                else
                    j = positions_mapping.synthetic[n,i]
                end

                # reactivities
                shape_reactivities_synth[i, n, c] = prof_df.HQ_profile[j]
                # reactivity standard errors
                shape_reactivities_synth_err[i, n, c] = prof_df.HQ_stderr[j]

                shape_raw_reactivities_synth[i,n,c] = prof_df.Reactivity_profile[j]
                shape_raw_reactivities_synth_err[i,n,c] = prof_df.Std_err[j]

                # mutation rates
                shape_M_synth[i,n,c] = prof_df.Modified_rate[j]
                shape_U_synth[i,n,c] = prof_df.Untreated_rate[j]
                shape_D_synth[i,n,c] = prof_df.Denatured_rate[j]

                # read depth
                shape_M_synth_depth[i,n,c] = prof_df.Modified_effective_depth[j]
                shape_U_synth_depth[i,n,c] = prof_df.Untreated_effective_depth[j]
                shape_D_synth_depth[i,n,c] = prof_df.Denatured_effective_depth[j]

                if "Norm_profile" ∈ names(prof_df) # some profiles do not have normalized reactivities because of few reads
                    shape_reactivities_indivnorm_synth[i,n,c] = prof_df.Norm_profile[j]
                    shape_reactivities_indivnorm_synth_err[i,n,c] = prof_df.Norm_stderr[j]
                end
            end
        end

        for n = 1:2 # PDB sequences
            aptamer_name = aptamer_names[306 + n]
            @assert aptamer_name ∈ ("SAMAP-PDB0", "SAMAP-PDB10")

            # _profile.txt files (no normalization!)
            prof_file = joinpath(shape_dir, cond, "$(cond)_$(aptamer_name)_profile.txt")

            index_in_aptamers_df = only(i for (i, name) = enumerate(aptamers_df.name) if occursin(aptamer_name, name))
            pdb_id = ("2GIS", "4KQY")[n]

            if isfile(prof_file)
                #@info "Loading $prof_file"
                prof_df = CSV.read(prof_file, DataFrame)
                @assert uppercase(join(prof_df.Sequence)) == replace(aptamers_df.sequence[index_in_aptamers_df], 'T' => 'U')
                for i = 1:108
                    if ismissing(positions_mapping.pdb[n,i])
                        continue
                    else
                        j = positions_mapping.pdb[n,i]
                    end

                    shape_reactivities_pdb[i, n, c] = prof_df.HQ_profile[j]
                    shape_reactivities_pdb_err[i, n, c] = prof_df.HQ_stderr[j]

                    shape_raw_reactivities_pdb[i, n, c] = prof_df.Reactivity_profile[j]
                    shape_raw_reactivities_pdb_err[i, n, c] = prof_df.Std_err[j]

                    shape_M_pdb[i,n,c] = prof_df.Modified_rate[j]
                    shape_U_pdb[i,n,c] = prof_df.Untreated_rate[j]
                    shape_D_pdb[i,n,c] = prof_df.Denatured_rate[j]

                    shape_M_pdb_depth[i,n,c] = prof_df.Modified_effective_depth[j]
                    shape_U_pdb_depth[i,n,c] = prof_df.Untreated_effective_depth[j]
                    shape_D_pdb_depth[i,n,c] = prof_df.Denatured_effective_depth[j]

                    if "Norm_profile" ∈ names(prof_df)
                        shape_reactivities_indivnorm_pdb[i,n,c] = prof_df.Norm_profile[j]
                        shape_reactivities_indivnorm_pdb_err[i,n,c] = prof_df.Norm_stderr[j]
                    end
                end
            else
                @warn "$prof_file not found!"
                push!(missing_files, prof_file)
            end
        end
    end

    # concat things
    # attention! I put natural first, then synthetic, then PDB, as in aptamers_df
    # NaN denotes missing values (due to alignment, or experimental noise)
    shape_reactivities = cat(shape_reactivities_natur, shape_reactivities_synth, shape_reactivities_pdb; dims=2)
    shape_reactivities_err = cat(shape_reactivities_natur_err, shape_reactivities_synth_err, shape_reactivities_pdb_err; dims=2)
    shape_raw_reactivities = cat(shape_raw_reactivities_natur, shape_raw_reactivities_synth, shape_raw_reactivities_pdb; dims=2)
    shape_raw_reactivities_err = cat(shape_raw_reactivities_natur_err, shape_raw_reactivities_synth_err, shape_raw_reactivities_pdb_err; dims=2)

    shape_reactivities_indivnorm = cat(shape_reactivities_indivnorm_natur, shape_reactivities_indivnorm_synth, shape_reactivities_indivnorm_pdb; dims=2)
    shape_reactivities_indivnorm_err = cat(shape_reactivities_indivnorm_natur_err, shape_reactivities_indivnorm_synth_err, shape_reactivities_indivnorm_pdb_err; dims=2)

    # mutation rates
    shape_M = cat(shape_M_natur, shape_M_synth, shape_M_pdb; dims=2)
    shape_U = cat(shape_U_natur, shape_U_synth, shape_U_pdb; dims=2)
    shape_D = cat(shape_D_natur, shape_D_synth, shape_D_pdb; dims=2)

    # read depth
    shape_M_depth = cat(shape_M_natur_depth, shape_M_synth_depth, shape_M_pdb_depth; dims=2)
    shape_U_depth = cat(shape_U_natur_depth, shape_U_synth_depth, shape_U_pdb_depth; dims=2)
    shape_D_depth = cat(shape_D_natur_depth, shape_D_synth_depth, shape_D_pdb_depth; dims=2)

    # standard error of the mutation rates (assuminng Poisson statistics)
    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    # consensus secondary structure in WUSS format of seed alignment
    wuss_full = stockholm_ss(Infernal.esl_afetch(Rfam.seed(), "RF00162").out)

    # sequence ids
    aptamer_ids = aptamers_df.id[indexin(aptamer_names, aptamers_df.name)]

    # fetch Rfam seed sequences
    RF00162_seed_afa = Infernal.esl_reformat(
        "AFA", Infernal.esl_afetch(Rfam.seed(), "RF00162").out; informat="STOCKHOLM"
    )
    RF00162_seed_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))
    RF00162_seed_seqs_full = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out))))
    RF00162_seed_match_cols = findall(≠('.'), wuss_full)
    RF00162_seed_seqs = [LongRNA{4}(seq[RF00162_seed_match_cols]) for seq in RF00162_seed_seqs_full]

    # fetch Rfam hits sequences -- already trimmed (no inserts) aligned fasta
    RF00162_hits_afa = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(), "RF00162").out,
        Rfam.fasta_file("RF00162");
        matchonly=true, outformat="AFA"
    )
    RF00162_hits_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out)))
    RF00162_hits_seqs = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))))

    RF00162_hits_afa_1410 = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out,
        Rfam.fasta_file("RF00162"; rfam_version="14.10");
        matchonly=true, outformat="afa"
     )
     RF00162_hits_dsc_1410 = FASTX.description.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out)))
     RF00162_hits_seqs_1410 = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out))))

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs_1410) .== 108)
    @assert length(RF00162_seed_match_cols) == 108

    ## Aligned sequences in indexed in the same order as the SHAPE data ... we will fill this array below.
    aligned_sequences = Vector{Union{Missing,String}}(undef, length(aptamer_names))

    # seed sequences (no missmatches, so no need for isnothing check)
    aligned_sequences[aptamer_origin .== "RF00162_seed70"] .= string.(
        RF00162_seed_seqs[indexin(aptamer_ids[aptamer_origin .== "RF00162_seed70"], RF00162_seed_ids)]
    )

    # hits (there are 5 missmatches, therefore we need the `isnothing` check)
    aligned_sequences[aptamer_origin .== "RF00162_full30"] .= [
        isnothing(i) ? missing : string(RF00162_hits_seqs[i]) for i = indexin(aptamer_ids[aptamer_origin .== "RF00162_full30"], RF00162_hits_ids)
    ]

    # file containing our synthetic sequences probed in 2022 (aligned!)
    __synth_df = probed_artificial_sequences_2022_df()

    # synthetic sequences
    for (i, j) in zip(1:100, indexin(__synth_df.Name, aptamer_names))
        if isnothing(j)
            aligned_sequences[j] = missing
        else
            @assert aptamer_origin[j] ∈ ("RF00162_syn_inf", "RF00162_syn_rbm")
            @assert aptamer_origin[j] == "RF00162_syn_" * __synth_df.origin[i]
            aligned_sequences[j] = __synth_df.sequence[i]
        end
    end

    # PDB sequences
    aligned_sequences[307] = string(RF00162_hits_seqs_1410[only(i for (i, dsc) = enumerate(RF00162_hits_dsc_1410) if occursin("2GIS", dsc))])
    aligned_sequences[308] = string(RF00162_hits_seqs_1410[only(i for (i, dsc) = enumerate(RF00162_hits_dsc_1410) if occursin("4KQY", dsc))])

    @assert all(ismissing.(aligned_sequences) .|| (length.(aligned_sequences) .== 108))
    @assert allunique(aptamer_names)
    @assert allunique(aptamer_ids)

    return (; shape_reactivities, shape_reactivities_err,
        shape_reactivities_indivnorm, shape_reactivities_indivnorm_err,
        shape_raw_reactivities, shape_raw_reactivities_err,
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        conditions, aptamer_names, aptamer_origin, aptamer_ids, aligned_sequences,
        missing_files
    )
end


function load_shapemapper_data_pierre_demux_20240801_with_pdb_repls_merged()
    shape_dir = load_shapemapper_data_pierre_demux_20240801_dir()
    conditions = filter(startswith("SAMAP"), readdir(shape_dir))
    @assert issorted(conditions)

    aptamers_df = probed_aptamers_table_20221027()
    synthetic_df = CSV.read(pierre20221107post(:synthetic), DataFrame)
    positions_mapping = shape_positions_alignment_2022()

    # Store reactivities in a 3-dimensional tensor,
    # shape_reactivities[i, n, c] = reactivity of site 'i', of sequence 'n', in condition 'c'
    # `missing` denotes sites without measurement (or deletions in the sequence)
    shape_reactivities_natur = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences
    shape_reactivities_pdb = fill(NaN, 108, 2, length(conditions)) # for PDB sequences
    shape_reactivities_natur_err = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth_err = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences
    shape_reactivities_pdb_err = fill(NaN, 108, 2, length(conditions)) # for PDB sequences

    # Reactivity profile after normalization, see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#reactivity-profile-calculation-and-normalization
    shape_reactivities_indivnorm_natur = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_pdb = fill(NaN, 108, 2, length(conditions))
    shape_reactivities_indivnorm_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth_err = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_pdb_err = fill(NaN, 108, 2, length(conditions))

    # Raw reactivities (without ShapeMapper quality checks)
    shape_raw_reactivities_natur = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_pdb = fill(NaN, 108, 2, length(conditions))
    shape_raw_reactivities_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth_err = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_pdb_err = fill(NaN, 108, 2, length(conditions))

    shape_M_natur = fill(NaN, 108, 206, length(conditions))
    shape_U_natur = fill(NaN, 108, 206, length(conditions))
    shape_D_natur = fill(NaN, 108, 206, length(conditions))
    shape_M_synth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth = fill(NaN, 108, 100, length(conditions))
    shape_M_pdb = fill(NaN, 108, 2, length(conditions))
    shape_U_pdb = fill(NaN, 108, 2, length(conditions))
    shape_D_pdb = fill(NaN, 108, 2, length(conditions))

    # these will contain effective read depths
    # see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#effective-read-depth
    shape_M_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_U_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_D_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_M_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_M_pdb_depth = fill(NaN, 108, 2, length(conditions))
    shape_U_pdb_depth = fill(NaN, 108, 2, length(conditions))
    shape_D_pdb_depth = fill(NaN, 108, 2, length(conditions))

    # aptamer names, consistent with the order in the shape_reactivities tensor
    aptamer_names = [
        ["APSAMN$n" for n in 1:206]; # natural
        ["APSAMS$n" for n in 1:100]; # synthetic
        ["SAMAP-PDB0", "SAMAP-PDB10"] # PDB
    ]
    # origin of aptamers (RF00162 natural, hits or seed, and the synthetic ones, RBM or Infernal)
    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
        ["RF00162_syn_" * synthetic_df.origin[n] for n in 1:100]; # synthetic (rbm or infernal)
        ["PDB", "PDB"] # from PDB
    ]
    @assert length(aptamer_names) == length(aptamer_origin) == 308

    missing_files = String[]

    for (c, cond) = enumerate(conditions)
        for n = 1:206 # natural sequences (from RF00162)
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

                    if "Norm_profile" ∈ names(prof_df)
                        shape_reactivities_indivnorm_natur[i,n,c] = prof_df.Norm_profile[j]
                        shape_reactivities_indivnorm_natur_err[i,n,c] = prof_df.Norm_stderr[j]
                    end
                end
            else
                @warn "$prof_file not found!"
                push!(missing_files, prof_file)
            end
        end

        for n = 1:100 # synthetic sequences
            prof_file = joinpath(shape_dir, cond, "$(cond)_APSAMS$(n)_profile.txt")
            @assert isfile(prof_file) "File $prof_file not found!"
            prof_df = CSV.read(prof_file, DataFrame)
            @assert uppercase(join(prof_df.Sequence)) == replace(aptamers_df.sequence[aptamers_df.source .== "RF00162_syn"][n], 'T' => 'U')
            for i in 1:108
                if ismissing(positions_mapping.synthetic[n,i])
                    continue
                else
                    j = positions_mapping.synthetic[n,i]
                end

                # reactivities
                shape_reactivities_synth[i, n, c] = prof_df.HQ_profile[j]
                # reactivity standard errors
                shape_reactivities_synth_err[i, n, c] = prof_df.HQ_stderr[j]

                shape_raw_reactivities_synth[i,n,c] = prof_df.Reactivity_profile[j]
                shape_raw_reactivities_synth_err[i,n,c] = prof_df.Std_err[j]

                # mutation rates
                shape_M_synth[i,n,c] = prof_df.Modified_rate[j]
                shape_U_synth[i,n,c] = prof_df.Untreated_rate[j]
                shape_D_synth[i,n,c] = prof_df.Denatured_rate[j]

                # read depth
                shape_M_synth_depth[i,n,c] = prof_df.Modified_effective_depth[j]
                shape_U_synth_depth[i,n,c] = prof_df.Untreated_effective_depth[j]
                shape_D_synth_depth[i,n,c] = prof_df.Denatured_effective_depth[j]

                if "Norm_profile" ∈ names(prof_df) # some profiles do not have normalized reactivities because of few reads
                    shape_reactivities_indivnorm_synth[i,n,c] = prof_df.Norm_profile[j]
                    shape_reactivities_indivnorm_synth_err[i,n,c] = prof_df.Norm_stderr[j]
                end
            end
        end

        for n = 1:2 # PDB sequences
            aptamer_name = aptamer_names[306 + n]
            @assert aptamer_name ∈ ("SAMAP-PDB0", "SAMAP-PDB10")

            # _profile.txt files (no normalization!)
            prof_file = joinpath(shape_dir, cond, "$(cond)_$(aptamer_name)_profile.txt")

            index_in_aptamers_df = only(i for (i, name) = enumerate(aptamers_df.name) if occursin(aptamer_name, name))
            pdb_id = ("2GIS", "4KQY")[n]

            if isfile(prof_file)
                #@info "Loading $prof_file"
                prof_df = CSV.read(prof_file, DataFrame)
                @assert uppercase(join(prof_df.Sequence)) == replace(aptamers_df.sequence[index_in_aptamers_df], 'T' => 'U')
                for i = 1:108
                    if ismissing(positions_mapping.pdb[n,i])
                        continue
                    else
                        j = positions_mapping.pdb[n,i]
                    end

                    shape_reactivities_pdb[i, n, c] = prof_df.HQ_profile[j]
                    shape_reactivities_pdb_err[i, n, c] = prof_df.HQ_stderr[j]

                    shape_raw_reactivities_pdb[i, n, c] = prof_df.Reactivity_profile[j]
                    shape_raw_reactivities_pdb_err[i, n, c] = prof_df.Std_err[j]

                    shape_M_pdb[i,n,c] = prof_df.Modified_rate[j]
                    shape_U_pdb[i,n,c] = prof_df.Untreated_rate[j]
                    shape_D_pdb[i,n,c] = prof_df.Denatured_rate[j]

                    shape_M_pdb_depth[i,n,c] = prof_df.Modified_effective_depth[j]
                    shape_U_pdb_depth[i,n,c] = prof_df.Untreated_effective_depth[j]
                    shape_D_pdb_depth[i,n,c] = prof_df.Denatured_effective_depth[j]

                    if "Norm_profile" ∈ names(prof_df)
                        shape_reactivities_indivnorm_pdb[i,n,c] = prof_df.Norm_profile[j]
                        shape_reactivities_indivnorm_pdb_err[i,n,c] = prof_df.Norm_stderr[j]
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
    shape_reactivities = cat(shape_reactivities_natur, shape_reactivities_synth, shape_reactivities_pdb; dims=2)
    shape_reactivities_err = cat(shape_reactivities_natur_err, shape_reactivities_synth_err, shape_reactivities_pdb_err; dims=2)
    shape_raw_reactivities = cat(shape_raw_reactivities_natur, shape_raw_reactivities_synth, shape_raw_reactivities_pdb; dims=2)
    shape_raw_reactivities_err = cat(shape_raw_reactivities_natur_err, shape_raw_reactivities_synth_err, shape_raw_reactivities_pdb_err; dims=2)

    shape_reactivities_indivnorm = cat(shape_reactivities_indivnorm_natur, shape_reactivities_indivnorm_synth, shape_reactivities_indivnorm_pdb; dims=2)
    shape_reactivities_indivnorm_err = cat(shape_reactivities_indivnorm_natur_err, shape_reactivities_indivnorm_synth_err, shape_reactivities_indivnorm_pdb_err; dims=2)

    # mutation rates
    shape_M = cat(shape_M_natur, shape_M_synth, shape_M_pdb; dims=2)
    shape_U = cat(shape_U_natur, shape_U_synth, shape_U_pdb; dims=2)
    shape_D = cat(shape_D_natur, shape_D_synth, shape_D_pdb; dims=2)

    # read depth
    shape_M_depth = cat(shape_M_natur_depth, shape_M_synth_depth, shape_M_pdb_depth; dims=2)
    shape_U_depth = cat(shape_U_natur_depth, shape_U_synth_depth, shape_U_pdb_depth; dims=2)
    shape_D_depth = cat(shape_D_natur_depth, shape_D_synth_depth, shape_D_pdb_depth; dims=2)

    # standard error of the mutation rates (assuminng Poisson statistics)
    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    # consensus secondary structure in WUSS format of seed alignment
    wuss_full = stockholm_ss(Infernal.esl_afetch(Rfam.seed(), "RF00162").out)

    # sequence ids
    aptamer_ids = aptamers_df.id[indexin(aptamer_names, aptamers_df.name)]

    # fetch Rfam seed sequences
    RF00162_seed_afa = Infernal.esl_reformat(
        "AFA", Infernal.esl_afetch(Rfam.seed(), "RF00162").out; informat="STOCKHOLM"
    )
    RF00162_seed_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))
    RF00162_seed_seqs_full = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_seed_afa.out))))
    RF00162_seed_match_cols = findall(≠('.'), wuss_full)
    RF00162_seed_seqs = [LongRNA{4}(seq[RF00162_seed_match_cols]) for seq in RF00162_seed_seqs_full]

    # fetch Rfam hits sequences -- already trimmed (no inserts) aligned fasta
    RF00162_hits_afa = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(), "RF00162").out,
        Rfam.fasta_file("RF00162");
        matchonly=true, outformat="AFA"
    )
    RF00162_hits_ids = FASTX.identifier.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out)))
    RF00162_hits_seqs = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))))

    RF00162_hits_afa_1410 = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out,
        Rfam.fasta_file("RF00162"; rfam_version="14.10");
        matchonly=true, outformat="afa"
     )
     RF00162_hits_dsc_1410 = FASTX.description.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out)))
     RF00162_hits_seqs_1410 = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out))))

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs_1410) .== 108)
    @assert length(RF00162_seed_match_cols) == 108

    ## Aligned sequences in indexed in the same order as the SHAPE data ... we will fill this array below.
    aligned_sequences = Vector{Union{Missing,String}}(undef, length(aptamer_names))

    # seed sequences (no missmatches, so no need for isnothing check)
    aligned_sequences[aptamer_origin .== "RF00162_seed70"] .= string.(
        RF00162_seed_seqs[indexin(aptamer_ids[aptamer_origin .== "RF00162_seed70"], RF00162_seed_ids)]
    )

    # hits (there are 5 missmatches, therefore we need the `isnothing` check)
    aligned_sequences[aptamer_origin .== "RF00162_full30"] .= [
        isnothing(i) ? missing : string(RF00162_hits_seqs[i]) for i in indexin(aptamer_ids[aptamer_origin .== "RF00162_full30"], RF00162_hits_ids)
    ]

    # file containing our synthetic sequences probed in 2022 (aligned!)
    __synth_df = probed_artificial_sequences_2022_df()

    # synthetic sequences
    for (i, j) in zip(1:100, indexin(__synth_df.Name, aptamer_names))
        if isnothing(j)
            aligned_sequences[j] = missing
        else
            @assert aptamer_origin[j] ∈ ("RF00162_syn_inf", "RF00162_syn_rbm")
            @assert aptamer_origin[j] == "RF00162_syn_" * __synth_df.origin[i]
            aligned_sequences[j] = __synth_df.sequence[i]
        end
    end

    # PDB sequences
    aligned_sequences[307] = string(RF00162_hits_seqs_1410[only(i for (i, dsc) = enumerate(RF00162_hits_dsc_1410) if occursin("2GIS", dsc))])
    aligned_sequences[308] = string(RF00162_hits_seqs_1410[only(i for (i, dsc) = enumerate(RF00162_hits_dsc_1410) if occursin("4KQY", dsc))])

    @assert all(ismissing.(aligned_sequences) .|| (length.(aligned_sequences) .== 108))
    @assert allunique(aptamer_names)
    @assert allunique(aptamer_ids)

    return (; shape_reactivities, shape_reactivities_err,
        shape_reactivities_indivnorm, shape_reactivities_indivnorm_err,
        shape_raw_reactivities, shape_raw_reactivities_err,
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        conditions, aptamer_names, aptamer_origin, aptamer_ids, aligned_sequences,
        missing_files
    )
end

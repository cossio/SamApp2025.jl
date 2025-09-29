function load_dms_data_20250303()
    dms_dir = artifact_dir_dms_data_2025_03_03()
    conditions = filter(startswith("SAMAP_DMS"), readdir(dms_dir))
    @assert issorted(conditions)

    # table of DMS tested sequences
    dms_df = load_dms_data_sequences_table_20250303_with_aligned_sequences()

    df_seqs_500 = designed_sequences_20230317()
    yann_sequences = map(FASTX.sequence, yann_sequences_20230322())
    df_2022 = probed_aptamers_table_20221027()

    positions_mapping_2022 = shape_positions_alignment_2022()

    # load previous SHAPE datasets for reference
    shape_data_500 = SamApp2025.load_shapemapper_data_500v2_20240315()
    shape_data_2023 = SamApp2025.load_shapemapper_data_pierre_demux_20230920()

    shape_data_2023_names = shape_data_2023.aptamer_names
    shape_data_500_names = ["APSAM-S2-" * lpad(n, 3, '0') for n = 0:499]
    @assert isempty(shape_data_2023_names ∩ shape_data_500_names)

    number_of_sequences = length(dms_df.name ∩ shape_data_2023_names) + length(dms_df.name ∩ shape_data_500_names)
    sequence_length = 108

    shape_M = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_U = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_D = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_M_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_U_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_D_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_reactivities = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_reactivities_err = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_raw_reactivities = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_raw_reactivities_err = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    aptamer_names = [name for name = dms_df.name if name ∈ shape_data_2023_names || name ∈ shape_data_500_names]
    aligned_sequence = [dms_df.aligned_sequence[i] for i = eachindex(dms_df.name) if isassigned(dms_df.aligned_sequence, i) && (dms_df.name[i] ∈ shape_data_2023_names || dms_df.name[i] ∈ shape_data_500_names)]
    @assert length(aptamer_names) == length(aligned_sequence) == number_of_sequences

    for (c, cond) = enumerate(conditions)
        for (n, name) = enumerate(aptamer_names)
            @assert name ∈ shape_data_2023_names || name ∈ shape_data_500_names

            profile_file = joinpath(dms_dir, cond, "$(cond)_$(name)_profile.txt")
            profile_df = CSV.read(profile_file, DataFrame)

            if name ∈ shape_data_2023_names
                index_in_2023 = only(findall(==(name), shape_data_2023_names))
                @assert uppercase(join(profile_df.Sequence)) == replace(df_2022.sequence[only(findall(df_2022.name .== name))], 'T' => 'U')

                for i = 1:108
                    @assert 1 ≤ index_in_2023 ≤ 306
                    if index_in_2023 ≤ 206 # natural
                        mapped_position = positions_mapping_2022.natural[index_in_2023, i]
                    else
                        mapped_position = positions_mapping_2022.synthetic[index_in_2023 - 206, i]
                    end

                    if ismissing(mapped_position)
                        continue
                    else
                        shape_M[i, n, c] = profile_df.Modified_rate[mapped_position]
                        shape_U[i, n, c] = profile_df.Untreated_rate[mapped_position]
                        shape_D[i, n, c] = profile_df.Denatured_rate[mapped_position]

                        shape_M_depth[i, n, c] = profile_df.Modified_effective_depth[mapped_position]
                        shape_U_depth[i, n, c] = profile_df.Untreated_effective_depth[mapped_position]
                        shape_D_depth[i, n, c] = profile_df.Denatured_effective_depth[mapped_position]

                        shape_reactivities[i, n, c] = profile_df.HQ_profile[mapped_position]
                        shape_reactivities_err[i, n, c] = profile_df.HQ_stderr[mapped_position]

                        shape_raw_reactivities[i, n, c] = profile_df.Reactivity_profile[mapped_position]
                        shape_raw_reactivities_err[i, n, c] = profile_df.Std_err[mapped_position]
                    end
                end
            elseif name ∈ shape_data_500_names
                index_in_500 = only(findall(==(name), shape_data_500_names))

                if index_in_500 ≤ 450 # our sequence
                    @assert uppercase(join(profile_df.Sequence)) == df_seqs_500.sequence[index_in_500]
                else # yann sequence
                    @assert uppercase(join(profile_df.Sequence)) == yann_sequences[index_in_500 - 450]
                end

                for i = 1:108
                    # for our sequences, consider only ungapped positions
                    index_in_500 ≤ 450 && df_seqs_500.aligned_sequence[index_in_500][i] == '-' && continue

                    j = i
                    if index_in_500 ≤ 450
                        algn = cumsum(collect(df_seqs_500.aligned_sequence[index_in_500]) .== '-')
                        j -= algn[i] # position in sequence without gaps
                    end

                    shape_M[i, n, c] = profile_df.Modified_rate[j]
                    shape_U[i, n, c] = profile_df.Untreated_rate[j]
                    shape_D[i, n, c] = profile_df.Denatured_rate[j]

                    shape_M_depth[i, n, c] = profile_df.Modified_effective_depth[j]
                    shape_U_depth[i, n, c] = profile_df.Untreated_effective_depth[j]
                    shape_D_depth[i, n, c] = profile_df.Denatured_effective_depth[j]

                    shape_reactivities[i, n, c] = profile_df.HQ_profile[j]
                    shape_reactivities_err[i, n, c] = profile_df.HQ_stderr[j]

                    shape_raw_reactivities[i, n, c] = profile_df.Reactivity_profile[j]
                    shape_raw_reactivities_err[i, n, c] = profile_df.Std_err[j]
                end
            else
                error("This should not happen")
            end
        end
    end

    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    return (;
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        shape_reactivities, shape_reactivities_err,
        shape_raw_reactivities, shape_raw_reactivities_err,
        conditions, aptamer_names, aligned_sequence
    )
end

function load_dms_data_sequences_20250303()
    file = joinpath(artifact_dir_dms_data_2025_03_03(), "2025-02-18-tagged-aptamer-sequences-SetA.fasta")
    sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(file)))
    return LongRNA{4}.([replace(seq, 'T' => 'U') for seq = sequences])
end

function load_dms_data_sequences_table_20250303()
    file = joinpath(artifact_dir_dms_data_2025_03_03(), "2025-02-18-tagged-aptamer-sequences-SetA.tsv")
    df = CSV.read(file, DataFrame)
    df.sequence = [replace(seq, 'T' => 'U') for seq = df.sequence]
    return df
end

function load_dms_data_sequences_table_20250303_with_aligned_sequences()
    df = load_dms_data_sequences_table_20250303()
    df.aligned_sequence = similar(df.sequence, Union{String, Missing})

    df_seqs_500 = designed_sequences_20230317()
    for (n, i) = enumerate(indexin(df.sequence, df_seqs_500.sequence))
        if isnothing(i)
            continue
        else
            df.aligned_sequence[n] = df_seqs_500.aligned_sequence[i]
        end
    end

    yann_seqs = map(FASTX.sequence, yann_sequences_20230322())
    for (n, i) = enumerate(indexin(df.sequence, yann_seqs))
        if isnothing(i)
            continue
        else
            df.aligned_sequence[n] = yann_seqs[i]
        end
    end

    shape_data_2023 = load_shapemapper_data_pierre_demux_20230920()
    for (n, i) = enumerate(indexin(df.name, shape_data_2023.aptamer_names))
        if isnothing(i)
            continue
        else
            df.aligned_sequence[n] = shape_data_2023.aligned_sequences[i]
        end
    end

     return df
end

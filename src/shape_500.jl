function load_shapemapper_data_500v2_20240315()
    shape_dir = artifact"2024-03-14-SAMAPS2-500Artificial-resequenced-condition"
    conditions = filter(startswith("SAMAPS2"), readdir(shape_dir))
    @assert issorted(conditions)

    our_design_df = designed_sequences_20230317()
    yann_sequences = map(FASTX.sequence, yann_sequences_20230322())

    shape_M = fill(NaN, 108, 500, length(conditions))
    shape_U = fill(NaN, 108, 500, length(conditions))
    shape_D = fill(NaN, 108, 500, length(conditions))

    shape_M_depth = fill(NaN, 108, 500, length(conditions))
    shape_U_depth = fill(NaN, 108, 500, length(conditions))
    shape_D_depth = fill(NaN, 108, 500, length(conditions))

    shape_reactivities = fill(NaN, 108, 500, length(conditions))
    shape_reactivities_err = fill(NaN, 108, 500, length(conditions))

    shape_raw_reactivities = fill(NaN, 108, 500, length(conditions))
    shape_raw_reactivities_err = fill(NaN, 108, 500, length(conditions))

    shape_reactivities_indivnorm = fill(NaN, 108, 500, length(conditions))
    shape_reactivities_indivnorm_err = fill(NaN, 108, 500, length(conditions))

    for (c, cond) in enumerate(conditions)
        for n = 1:450 # our sequences
            profile_file = joinpath(shape_dir, cond, "$(cond)_APSAM-S2-$(lpad(n - 1, 3, '0'))_profile.txt")
            profile_df = CSV.read(profile_file, DataFrame)

            if n ≤ 450
                # our sequence
                @assert uppercase(join(profile_df.Sequence)) == our_design_df.sequence[n]
            else
                # n > 450, Yann's sequence
                @assert uppercase(join(profile_df.Sequence)) == yann_sequences[n - 450]
            end

            for i = 1:108
                # for our sequences, consider only ungapped positions
                n ≤ 450 && our_design_df.aligned_sequence[n][i] == '-' && continue

                j = i

                if n ≤ 450
                    algn = cumsum(collect(our_design_df.aligned_sequence[n]) .== '-')
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

                # some profiles do not have normalized reactivities because of few reads
                if "Norm_profile" ∈ names(profile_df)
                    shape_reactivities_indivnorm[i, n, c] = profile_df.Norm_profile[j]
                    shape_reactivities_indivnorm_err[i, n, c] = profile_df.Norm_stderr[j]
                end
            end
        end
    end

    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    aptamer_origin = [our_design_df.origin; fill("Infrared", 50)]
    aptamer_criteria = [our_design_df.criteria; fill("-", 50)]
    aligned_sequences = [our_design_df.aligned_sequence; yann_sequences]

    return (; shape_reactivities, shape_reactivities_err,
        shape_reactivities_indivnorm,
        shape_reactivities_indivnorm_err,
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        shape_raw_reactivities, shape_raw_reactivities_err,
        conditions, aptamer_origin, aptamer_criteria,
        aligned_sequences = LongRNA{4}.(aligned_sequences)
    )
end

#= File containing our designed sequences for second run of experiments (with the 500 aptamers) =#
function designed_sequences_20230317()
    CSV.read(joinpath(artifact"Designed_Sequences_20230317", "all_designed_sequences_v3.csv"), DataFrame)
end

function yann_sequences_20230322()
    fasta = joinpath(artifact"yann_sequences_20230322", "yann_sequences", "rational-designs-aptamer.fa")
    return collect(FASTX.FASTA.Reader(open(fasta)))
end

const hallmark_sites_20230507 = [
    10, 11, # P1
    25, 26, 27, 28, 77, 79, # Pseudoknot
    34, 35, 36, 37, # Kink-turn
    73, 74, 75, # A-minor
    46, 47, # SAM contact (bulge)
    76, 100, # Base-triple
    101, 102, 103, 104, 105 # SAM contact (P1)
];



# a new version of the select_conditions function, adapted to the new namedtuple format
function select_conditions_20231002(shape_data, conditions)
    c = indexin(conditions, shape_data.conditions)
    return (;
        shape_data...,
        shape_reactivities = shape_data.shape_reactivities[:,:,c],
        shape_reactivities_err = shape_data.shape_reactivities_err[:,:,c],
        shape_raw_reactivities = shape_data.shape_raw_reactivities[:,:,c],
        shape_raw_reactivities_err = shape_data.shape_raw_reactivities_err[:,:,c],
        shape_reactivities_indivnorm = shape_data.shape_reactivities_indivnorm[:,:,c],
        shape_reactivities_indivnorm_err = shape_data.shape_reactivities_indivnorm_err[:,:,c],
        shape_M = shape_data.shape_M[:,:,c],
        shape_U = shape_data.shape_U[:,:,c],
        shape_D = shape_data.shape_D[:,:,c],
        shape_M_depth = shape_data.shape_M_depth[:,:,c],
        shape_U_depth = shape_data.shape_U_depth[:,:,c],
        shape_D_depth = shape_data.shape_D_depth[:,:,c],
        shape_M_stderr = shape_data.shape_M_stderr[:,:,c],
        shape_U_stderr = shape_data.shape_U_stderr[:,:,c],
        shape_D_stderr = shape_data.shape_D_stderr[:,:,c],
        # shape_gamma_theta_M = shape_data.shape_gamma_theta_M[:,:,c],
        # shape_gamma_theta_U = shape_data.shape_gamma_theta_U[:,:,c],
        # shape_gamma_theta_D = shape_data.shape_gamma_theta_D[:,:,c],
        # shape_gamma_alpha_M = shape_data.shape_gamma_alpha_M[:,:,c],
        # shape_gamma_alpha_U = shape_data.shape_gamma_alpha_U[:,:,c],
        # shape_gamma_alpha_D = shape_data.shape_gamma_alpha_D[:,:,c],
        conditions
    )
end

function shapemapper_data_pierre_demux_20230920_dir(; demux::Bool=true)
    @assert demux
    if demux
        shape_dir = artifact"SAMAP-demux-repl0-and-repl4+5-2023-09"
    else
        #shape_dir = artifact"SAMAP-nodemux-repl0-and-repl4+5-2023-09"
    end
    @assert isdir(shape_dir)
    return shape_dir
end

#=
Loads the SHAPE data produced by SHAPE-Mapper2 pipeline.
This part does not do the statistical analysis yet.
Data demuxed by Pierre, for replicates 4+5
Also available at "pCloud Drive/data/2023/SAM-Riboswitch/20230920-Pierre_demux"
(Here downloaded automatically as Julia Artifacts)
=#
function load_shapemapper_data_pierre_demux_20230920(; demux::Bool=true)
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
    shape_reactivities_natur_err = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth_err = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences

    # Reactivity profile after normalization, see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#reactivity-profile-calculation-and-normalization
    shape_reactivities_indivnorm_natur = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth_err = fill(NaN, 108, 100, length(conditions))

    # Raw reactivities (without ShapeMapper quality checks)
    shape_raw_reactivities_natur = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth_err = fill(NaN, 108, 100, length(conditions))

    shape_M_natur = fill(NaN, 108, 206, length(conditions))
    shape_U_natur = fill(NaN, 108, 206, length(conditions))
    shape_D_natur = fill(NaN, 108, 206, length(conditions))
    shape_M_synth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth = fill(NaN, 108, 100, length(conditions))

    # these will contain effective read depths
    # see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#effective-read-depth
    shape_M_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_U_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_D_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_M_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth_depth = fill(NaN, 108, 100, length(conditions))

    # aptamer names, consistent with the order in the shape_reactivities tensor
    aptamer_names = [
        ["APSAMN$n" for n in 1:206]; # natural
        ["APSAMS$n" for n in 1:100]  # synthetic
    ]
    # origin of aptamers (RF00162 natural, hits or seed, and the synthetic ones, RBM or Infernal)
    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
        ["RF00162_syn_" * synthetic_df.origin[n] for n in 1:100] # synthetic (rbm or infernal)
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

        for n in 1:100 # synthetic sequences
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
    end

    # concat things
    # attention! I put natural first, then synthetic, as in aptamers_df
    # NaN denotes missing values (due to alignment, or experimental noise)
    shape_reactivities = cat(shape_reactivities_natur, shape_reactivities_synth; dims=2)
    shape_reactivities_err = cat(shape_reactivities_natur_err, shape_reactivities_synth_err; dims=2)
    shape_raw_reactivities = cat(shape_raw_reactivities_natur, shape_raw_reactivities_synth; dims=2)
    shape_raw_reactivities_err = cat(shape_raw_reactivities_natur_err, shape_raw_reactivities_synth_err; dims=2)

    shape_reactivities_indivnorm = cat(shape_reactivities_indivnorm_natur, shape_reactivities_indivnorm_synth; dims=2)
    shape_reactivities_indivnorm_err = cat(shape_reactivities_indivnorm_natur_err, shape_reactivities_indivnorm_synth_err; dims=2)

    # mutation rates
    shape_M = cat(shape_M_natur, shape_M_synth; dims=2)
    shape_U = cat(shape_U_natur, shape_U_synth; dims=2)
    shape_D = cat(shape_D_natur, shape_D_synth; dims=2)

    # read depth
    shape_M_depth = cat(shape_M_natur_depth, shape_M_synth_depth; dims=2)
    shape_U_depth = cat(shape_U_natur_depth, shape_U_synth_depth; dims=2)
    shape_D_depth = cat(shape_D_natur_depth, shape_D_synth_depth; dims=2)

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

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
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


function aligned_sequences_20230920()
    synthetic_df = CSV.read(pierre20221107post(:synthetic), DataFrame)
    aptamers_df = probed_aptamers_table_20221027()

    aptamer_names = [
        ["APSAMN$n" for n in 1:206]; # natural
        ["APSAMS$n" for n in 1:100]  # synthetic
    ]

    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
        ["RF00162_syn_" * synthetic_df.origin[n] for n in 1:100] # synthetic (rbm or infernal)
    ]
    @assert length(aptamer_names) == length(aptamer_origin)

    wuss_full = stockholm_ss(Infernal.esl_afetch(Rfam.seed(), "RF00162").out)

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

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
    @assert length(RF00162_seed_match_cols) == 108

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

    @assert all(ismissing.(aligned_sequences) .|| (length.(aligned_sequences) .== 108))
    @assert allunique(aptamer_names)
    @assert allunique(aptamer_ids)

    return aligned_sequences
end



# Compute Protection scores (log-odd ratios)
function shape_basepair_log_odds_v4(;
    nsamples::Int=1000,
    p_thresh::Real=1e-3,
    shape_data, # aligned shape data

    # Reactivities assumed to be compatible with base-paired and unpaired states. We use
    # these to estimate the kernel density estimators for the probability of being
    # base-paired vs. unpaired. These should be normalized before!
    paired_reactivities::AbstractArray,
    unpaired_reactivities::AbstractArray,
    all_reactivities::AbstractArray,

    # normalization of reactivities (by default = 1)
    normalization::AbstractArray = one.(shape_data.shape_reactivities),

    only_hq_profile::Bool=false
)
    # now we estimate probabilities of being base-paired vs. unpaired based on these reactivities
    L, N, C = size(shape_data.shape_reactivities)

    # parameters of the posterior Gamma distribution (Gamma is the conjugate to the Poisson likelihood)
    # posterior Gamma distribution parameters (θ, α)
    shape_gamma_theta_M = (shape_data.shape_M_stderr.^2) ./ shape_data.shape_M
    shape_gamma_theta_U = (shape_data.shape_U_stderr.^2) ./ shape_data.shape_U
    shape_gamma_theta_D = (shape_data.shape_D_stderr.^2) ./ shape_data.shape_D

    shape_gamma_alpha_M = (shape_data.shape_M ./ shape_data.shape_M_stderr).^2
    shape_gamma_alpha_U = (shape_data.shape_U ./ shape_data.shape_U_stderr).^2
    shape_gamma_alpha_D = (shape_data.shape_D ./ shape_data.shape_D_stderr).^2

    # Kernel density estimator of PDFs of base-paired and unpaired reactivities
    bp_reactivities_kde = kde_lscv(filter(isfinite, vec(paired_reactivities))) # basepairs
    np_reactivities_kde = kde_lscv(filter(isfinite, vec(unpaired_reactivities))) # unpaired
    all_reactivities_kde = kde_lscv(filter(isfinite, vec(all_reactivities))) # ALL reactivities (base-paired + unpaired)

    # optimization to allow faster interpolation evaluation
    bp_reactivities_ikde = InterpKDE(bp_reactivities_kde)
    np_reactivities_ikde = InterpKDE(np_reactivities_kde)
    all_reactivities_ikde = InterpKDE(all_reactivities_kde)

    # we will fill these arrays below. These are the log-likelihood ratios of reactivities
    # conditioned on the site being base-paired (bp) or unpaired (np).
    # '0' refers to the "null" distribution of all reactivities together (all_reactivities_kde above)
    shape_log_odds_bp = fill(NaN, L, N, C) # log( P(R|bp) / P(R|0) )
    shape_log_odds_np = fill(NaN, L, N, C) # log( P(R|np) / P(R|0) )

    for c = 1:C, n = 1:N, i = 1:L
        α_M = shape_gamma_alpha_M[i,n,c]
        α_U = shape_gamma_alpha_U[i,n,c]
        α_D = shape_gamma_alpha_D[i,n,c]

        θ_M = shape_gamma_theta_M[i,n,c]
        θ_U = shape_gamma_theta_U[i,n,c]
        θ_D = shape_gamma_theta_D[i,n,c]

        # if something is NaN, just skip this site
        if isnan(α_M) || isnan(α_U) || isnan(α_D)
            continue
        elseif isnan(θ_M) || isnan(θ_U) || isnan(θ_D)
            continue
        end

        # discard if ShapeMapper rules think this site has poor statistics
        if only_hq_profile && isnan(shape_data.shape_reactivities[i,n,c])
            continue
        end

        # resample mutation rates
        M = rand(Gamma(α_M, θ_M), nsamples)
        U = rand(Gamma(α_U, θ_U), nsamples)
        D = rand(Gamma(α_D, θ_D), nsamples)

        # resampled reactivities
        R = (M - U) ./ D
        # normalize
        R = R / normalization[i,n,c]

        P_bp = pdf.(Ref(bp_reactivities_ikde), R)
        P_np = pdf.(Ref(np_reactivities_ikde), R)
        #P_0 = (P_bp + P_np) / 2 # this is the difference with respect to v2
        P_0 = pdf.(Ref(all_reactivities_ikde), R)

        # Apply a probability threshold (p_thresh) to remove outliers
        (mean(P_bp) > p_thresh) && (mean(P_np) > p_thresh) && (mean(P_0) > p_thresh) || continue
        # (mean(P_0) < p_thresh || mean(P_bp) < p_thresh && mean(P_np) < p_thresh) && continue

        _nz = findall(P_0 .> p_thresh) # retain only things that are not too far from the support of the histogram
        #@show mean(P_bp[_nz] ./ P_0[_nz]) mean(P_np[_nz] ./ P_0[_nz])

        shape_log_odds_bp[i,n,c] = log(mean(P_bp[_nz] ./ P_0[_nz]))
        shape_log_odds_np[i,n,c] = log(mean(P_np[_nz] ./ P_0[_nz]))
        @assert isfinite(shape_log_odds_bp[i,n,c])
        @assert isfinite(shape_log_odds_np[i,n,c])
    end

    # log( P(R|bp) / P(R|np) )
    shape_log_odds = shape_log_odds_bp - shape_log_odds_np

    return (; shape_data..., shape_log_odds, shape_log_odds_bp, shape_log_odds_np,
        shape_gamma_theta_M, shape_gamma_alpha_M,
        shape_gamma_theta_U, shape_gamma_alpha_U,
        shape_gamma_theta_D, shape_gamma_alpha_D
    )
end


function load_shapemapper_data_pierre_demux_20240801_dir()
    shape_dir = artifact"SAMAP_ALL-REP-MERGED-2023-10-27"
    @assert isdir(shape_dir)
    return shape_dir
end


#=
Merged replicates. Data demuxed by Pierre, and then processed by ShapeMapper.
=#
function load_shapemapper_data_pierre_demux_20231027_repls_merged()
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
    shape_reactivities_natur_err = fill(NaN, 108, 206, length(conditions)) # for natural sequences
    shape_reactivities_synth_err = fill(NaN, 108, 100, length(conditions)) # for synthetic sequences

    # Reactivity profile after normalization, see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#reactivity-profile-calculation-and-normalization
    shape_reactivities_indivnorm_natur = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth = fill(NaN, 108, 100, length(conditions))
    shape_reactivities_indivnorm_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_reactivities_indivnorm_synth_err = fill(NaN, 108, 100, length(conditions))

    # Raw reactivities (without ShapeMapper quality checks)
    shape_raw_reactivities_natur = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth = fill(NaN, 108, 100, length(conditions))
    shape_raw_reactivities_natur_err = fill(NaN, 108, 206, length(conditions))
    shape_raw_reactivities_synth_err = fill(NaN, 108, 100, length(conditions))

    shape_M_natur = fill(NaN, 108, 206, length(conditions))
    shape_U_natur = fill(NaN, 108, 206, length(conditions))
    shape_D_natur = fill(NaN, 108, 206, length(conditions))
    shape_M_synth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth = fill(NaN, 108, 100, length(conditions))

    # these will contain effective read depths
    # see https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/analysis_steps.md#effective-read-depth
    shape_M_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_U_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_D_natur_depth = fill(NaN, 108, 206, length(conditions))
    shape_M_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_U_synth_depth = fill(NaN, 108, 100, length(conditions))
    shape_D_synth_depth = fill(NaN, 108, 100, length(conditions))

    # aptamer names, consistent with the order in the shape_reactivities tensor
    aptamer_names = [
        ["APSAMN$n" for n in 1:206]; # natural
        ["APSAMS$n" for n in 1:100]  # synthetic
    ]
    # origin of aptamers (RF00162 natural, hits or seed, and the synthetic ones, RBM or Infernal)
    aptamer_origin = [
        ["RF00162_full30" for _ in 1:55];   # all hits of the alignment
        ["RF00162_seed70" for _ in 56:206]; # seed sequences in the alignment
        ["RF00162_syn_" * synthetic_df.origin[n] for n in 1:100] # synthetic (rbm or infernal)
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

        for n in 1:100 # synthetic sequences
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
    end

    # concat things
    # attention! I put natural first, then synthetic, as in aptamers_df
    # NaN denotes missing values (due to alignment, or experimental noise)
    shape_reactivities = cat(shape_reactivities_natur, shape_reactivities_synth; dims=2)
    shape_reactivities_err = cat(shape_reactivities_natur_err, shape_reactivities_synth_err; dims=2)
    shape_raw_reactivities = cat(shape_raw_reactivities_natur, shape_raw_reactivities_synth; dims=2)
    shape_raw_reactivities_err = cat(shape_raw_reactivities_natur_err, shape_raw_reactivities_synth_err; dims=2)

    shape_reactivities_indivnorm = cat(shape_reactivities_indivnorm_natur, shape_reactivities_indivnorm_synth; dims=2)
    shape_reactivities_indivnorm_err = cat(shape_reactivities_indivnorm_natur_err, shape_reactivities_indivnorm_synth_err; dims=2)

    # mutation rates
    shape_M = cat(shape_M_natur, shape_M_synth; dims=2)
    shape_U = cat(shape_U_natur, shape_U_synth; dims=2)
    shape_D = cat(shape_D_natur, shape_D_synth; dims=2)

    # read depth
    shape_M_depth = cat(shape_M_natur_depth, shape_M_synth_depth; dims=2)
    shape_U_depth = cat(shape_U_natur_depth, shape_U_synth_depth; dims=2)
    shape_D_depth = cat(shape_D_natur_depth, shape_D_synth_depth; dims=2)

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

    @assert all(length.(RF00162_seed_seqs) .== 108)
    @assert all(length.(RF00162_hits_seqs) .== 108)
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

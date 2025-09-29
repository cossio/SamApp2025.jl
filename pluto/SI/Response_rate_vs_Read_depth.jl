### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 953d7355-78c1-4309-9f5f-d0c02abffc51
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ bb8a6eeb-32a3-46dc-982a-d826253291e3
using BioSequences: LongRNA

# ╔═╡ 54d6f1c4-bb92-457c-94ec-e67dd06fd5d0
using DataFrames: DataFrame

# ╔═╡ 80e407e4-f584-4253-b685-28e92b5c9119
using Distributions: Gamma, Poisson

# ╔═╡ 29da1d84-e6bc-4890-9a02-9c49e6505681
using LinearAlgebra: Diagonal, eigen

# ╔═╡ 9bee752e-9272-46af-8bc0-8c2a2a1807e2
using Makie: @L_str

# ╔═╡ 53a19766-1d66-43ee-b6df-3ef581746a78
using NaNStatistics: nansum, nanmean

# ╔═╡ 8c65263e-1e8e-40d8-b824-915b98a5ce57
using Random: bitrand

# ╔═╡ 75f00c6e-2859-4b9e-8e6e-cb102f2fd071
using Statistics: cor, mean

# ╔═╡ 2486a3e8-6e80-4d84-a1b3-e8952f3f7faa
using StatsBase: countmap, corspearman

# ╔═╡ 1c6856dd-20cb-449f-abcb-545316b28ab5
md"# Imports"

# ╔═╡ 949fc97f-e30d-4b6c-aa9a-1ec7f1275df6
import Makie, CairoMakie, Logomaker, PlutoUI, KernelDensity, StatsBase, HDF5, CSV, FASTX, Infernal, Rfam, SamApp2025

# ╔═╡ 2e51c390-d4ee-4362-b7cb-409d7cdc8df8
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 571ede6a-ad6b-4710-a17d-9c3c7951d719
PlutoUI.TableOfContents()

# ╔═╡ a768ec77-2f9f-4e18-b488-b0cb756a4b00
md"# General data"

# ╔═╡ 735efa70-2f53-45e0-a301-2ea716ea242a
(; bps, nps, pks) = SamApp2025.RF00162_sites_paired()

# ╔═╡ aea6c6c9-faae-4c7e-bf35-b1faace52177
ENV["JULIA_RFAM_DIR"]

# ╔═╡ 84389a19-1f13-4162-a51b-19dcf16ea0c1
# Hallmark sites (without 10, 11)
_sites = SamApp2025.hallmark_sites_20230507[3:end]

# ╔═╡ ee44f513-ffa0-4183-bcf9-fdbcb73dba69
# signifcance threshold for protection scores
_thresh = log(5)

# ╔═╡ 781bdb61-4311-4e5c-95d4-4deb751c7e6e
md"# Load data"

# ╔═╡ aae30131-4654-4f98-83f1-9d0818189785
md"## Load data (Repl.0)"

# ╔═╡ 6781ca9d-85fc-46c9-8059-e62e13edabd2
# load SHAPE data
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ ed515a86-9939-4123-9237-7747817b03d0
# split rep0 from rep4+5
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ b222c039-9c93-4bcb-a4cb-bc4219b11f8f
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 392b2a08-597f-47e3-aeba-b45badcd43f5
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ a44fce13-f88e-456b-b870-983503331d16
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 341b5e17-056c-478d-8c3e-264bbf4bfeb0
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 1294fa8e-d76b-406c-8d07-21e3aac0c59c
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 6e36ec71-a706-4dab-be55-f2082f4a21c3
inf_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 0afc5cf7-852d-4d05-b3ad-7477bf86c7b9
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 7f93a86c-17e4-4472-a9db-c6b43fce0dca
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 4ff40840-74cb-4ac4-80b7-982c57478969
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ 600a0427-ad61-4317-aed9-41ba609c40c0
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 57444eac-7b7b-47c5-8e85-1ca01b4230ce
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 32abd57f-71b4-4188-9f28-0497c5ab7764
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 6508d862-75c5-4076-9aec-29043a77493e
shape_stats_rep0 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ daa2c7e2-2e6a-4a37-b32d-8569f21800bf
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3));

# ╔═╡ 5f1bb0e2-06e9-4476-8754-6175ed708b2e
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3));

# ╔═╡ a7484a8b-6feb-4ccc-987b-d87f800970f2
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 430ded66-dc00-411a-9db7-130f2d6555bf
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ e1b3f8ef-9b44-4bec-8610-01145774be1a
aptamer_rbm_energies_rep0 = [
    ismissing(seq) ? missing :
    RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ b0b44ab1-bd99-4239-b93f-3dea79d0c053
md"## Load data (500 seqs)"

# ╔═╡ a94baef6-c552-4b33-9139-0dc61dcf5aa9
shape_data_500 = SamApp2025.load_shapemapper_data_500v2_20240315();

# ╔═╡ 3cfc5665-936a-4621-8aa4-820c9a3e2c78
conds_sam_500, conds_mg_500, conds_30C_500 = [1,2], [4], [6];

# ╔═╡ 2e606abc-4ae4-4798-88ff-48153e41038b
bps_reactivities_500 = shape_data_500.shape_reactivities[bps, :, conds_sam_500];

# ╔═╡ d1b68ed8-f33e-4a02-b84f-03130ab96e58
nps_reactivities_500 = shape_data_500.shape_reactivities[nps, :, conds_sam_500];

# ╔═╡ 52a1f8c4-013a-4215-9b12-98f2dc23b2f3
all_reactivities_500 = shape_data_500.shape_reactivities[:, :, conds_sam_500];

# ╔═╡ f2cd1fd6-3b73-4615-9d3f-078dcdf50911
shape_stats_500 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = bps_reactivities_500,
    unpaired_reactivities = nps_reactivities_500,
    all_reactivities = all_reactivities_500,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ d51d3406-0e89-4a36-9591-199e3675d10f
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :,  conds_mg_500]; dim=(1,3));

# ╔═╡ b9badf2d-b58e-4b12-919c-feb79f911926
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ d192c0af-8551-4cba-b3c3-e8ea98c1b2a3
_responds_sam_yes_500 = (x_mg_500 .< -_thresh) .& (x_sam_500 .> +_thresh);

# ╔═╡ cb4d2651-8471-4fe6-ae96-413306f9ecc8
_responds_sam_nop_500 = (x_mg_500 .> +_thresh) .| (x_sam_500 .< -_thresh);

# ╔═╡ 1639c8dc-e96d-4225-bcda-b838ad391d8b
aptamer_rbm_energies_500 = [
    ismissing(seq) ? missing :
    RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq)))
    for seq in shape_data_500.aligned_sequences
];

# ╔═╡ 67ee804e-8877-4fc1-b2b5-a213bce545cd
md"# Make table"

# ╔═╡ 7e76da0f-91ca-4c5b-bf1d-ca10a2a6c5a8
_responsive_sam_rep0 = ifelse.(_responds_sam_yes_rep0, "Responsive", ifelse.(_responds_sam_nop_rep0, "Non-responsive", "Inconclusive"));

# ╔═╡ bfe80af1-a0ad-4702-9571-869aec104c28
_responsive_sam_500 = ifelse.(_responds_sam_yes_500, "Responsive", ifelse.(_responds_sam_nop_500, "Non-responsive", "Inconclusive"));

# ╔═╡ b5b6abab-40b0-430e-a012-59af58325a4d
ss = SamApp2025.RF00162_sites_annotated_secondary_structure();

# ╔═╡ 400df1dd-979a-45f4-9884-a6d9c012bd32
seq_groups_dfs = SamApp2025.artifact_load_sequencing_groups_2024_11_27();

# ╔═╡ dd49243c-3bbb-40fa-ac72-308c4ddebb57
group_names = sort(collect(keys(seq_groups_dfs)))

# ╔═╡ f6bb72e2-267f-45e6-b2f9-905bda307b06
begin
	df = DataFrame(;
		aptamer_names = [shape_data_rep0.aptamer_names; ["APSAM-S2-" * lpad(n - 1, 3, "0") for n = 1:500][shape_data_500.aptamer_origin .!= "Infrared"]],
	    aligned_sequences = [[ismissing(seq) ? missing : string(seq) for seq = shape_data_rep0.aligned_sequences]; string.(shape_data_500.aligned_sequences[shape_data_500.aptamer_origin .!= "Infrared"])],
	    aptamer_origin = [shape_data_rep0.aptamer_origin; shape_data_500.aptamer_origin[shape_data_500.aptamer_origin .!= "Infrared"]],
		experiment = [fill("Experiment_1", length(shape_data_rep0.aligned_sequences)); fill("Experiment_2", length(shape_data_500.aligned_sequences))[shape_data_500.aptamer_origin .!= "Infrared"]],
	    responsive = [_responsive_sam_rep0; _responsive_sam_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		RBM_score = [-aptamer_rbm_energies_rep0; -aptamer_rbm_energies_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		Protect_Score_Hallmark_Mg = [x_mg_rep0; x_mg_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		Protect_Score_Hallmark_SAM = [x_sam_rep0; x_sam_500[shape_data_500.aptamer_origin .!= "Infrared"]]
	)

	df.sequencing_group = [group_names[only(findall([n ∈ seq_groups_dfs[k].name for k = group_names]))] for n = df.aptamer_names]

	@assert count(df.experiment .== "Experiment_1") == length(shape_data_045.aligned_sequences) == 306
	for (seq1, seq2) = zip(df.aligned_sequences[df.experiment .== "Experiment_1"], shape_data_045.aligned_sequences)
		@assert ismissing(seq1) && ismissing(seq2) || seq1 == seq2
	end
	@assert df.experiment == [fill("Experiment_1", 306); fill("Experiment_2", 450)]

	df.read_depth_M = [
		[shape_data_rep0.shape_M_depth[:, n, :] for n = axes(shape_data_rep0.shape_M_depth, 2)];
		[shape_data_500.shape_M_depth[:, n, :] for n = 1:450]
	]
end

# ╔═╡ afa683f5-4bc3-4bfd-8e73-ba3e4c911368
let fig = Makie.Figure()
	primer_colors = (:orange, :purple, :red, :blue, :teal)
	#markers = (:utriangle, :xcross, :dtriangle, :circle, :cross)
	markers = (:circle, :circle, :circle, :circle, :circle)

	ax = Makie.Axis(fig[1,1]; width=230, height=230, xlabel="Read depth (M)", ylabel="Inconclusives / Total")
	for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
		primer = parse(Int, last(gr))

		read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[df.sequencing_group .== gr]))
		inconclusive_rate = mean(df.responsive[df.sequencing_group .== gr] .== "Inconclusive")
		if gr_n ≤ 5
			Makie.scatter!(ax, read_depth, inconclusive_rate; color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
		elseif gr_n < 9 # skip last one which has a single sequence
			Makie.scatter!(ax, read_depth, inconclusive_rate; color=primer_colors[primer], marker=markers[primer], markersize=15)
		end
	end
	Makie.xlims!(ax, 0, 3.3e4)
	Makie.ylims!(ax, -0.05, 0.65)

	fig[1,2] = Makie.Legend(fig, ax, "Primers", framevisible = false)

	for (col, (origin, title)) = enumerate(zip([["RF00162_full30", "RF00162_seed70"], ["RF00162_syn_rbm", "rbm"], ["RF00162_syn_inf", "infernal"]], ["Naturals", "RBM", "CM"]))
		ax = Makie.Axis(fig[1,3][1,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. / Total", title, xticks=[1e4, 3e4])
		for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
			_flag = (df.sequencing_group .== gr) .&& (df.aptamer_origin .∈ Ref(origin))
			primer = parse(Int, last(gr))
			if any(_flag)
				read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[_flag]))
				resp_rate = mean(df.responsive[_flag] .== "Responsive")
				if gr_n ≤ 5
					Makie.scatter!(ax, read_depth, resp_rate; color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
				elseif gr_n < 9 # skip last one which has a single sequence
					Makie.scatter!(ax, read_depth, resp_rate; color=primer_colors[primer], marker=markers[primer], markersize=15)
				end
			end
		end
		Makie.xlims!(ax, 0, 3.3e4)
		Makie.ylims!(ax, -0.05, 0.75)

		ax = Makie.Axis(fig[1,3][2,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. rate", xticks=[1e4, 3e4])
		for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
			_flag = (df.sequencing_group .== gr) .&& (df.aptamer_origin .∈ Ref(origin))
			if any(_flag)
				primer = parse(Int, last(gr))
				read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[_flag]))
				resp_rate = mean(df.responsive[_flag] .== "Responsive")
				inconclusive_rate = mean(df.responsive[df.sequencing_group .== gr] .== "Inconclusive")
				if gr_n ≤ 5
					Makie.scatter!(ax, read_depth, resp_rate / (1 - inconclusive_rate); color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
				elseif gr_n < 9 # skip last one which has a single sequence
					Makie.scatter!(ax, read_depth, resp_rate / (1 - inconclusive_rate); color=primer_colors[primer], marker=markers[primer], markersize=15)
				end
			end
		end
		Makie.xlims!(ax, 0, 3.3e4)
		Makie.ylims!(ax, -0.05, 0.75)

	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═1c6856dd-20cb-449f-abcb-545316b28ab5
# ╠═953d7355-78c1-4309-9f5f-d0c02abffc51
# ╠═949fc97f-e30d-4b6c-aa9a-1ec7f1275df6
# ╠═2e51c390-d4ee-4362-b7cb-409d7cdc8df8
# ╠═bb8a6eeb-32a3-46dc-982a-d826253291e3
# ╠═54d6f1c4-bb92-457c-94ec-e67dd06fd5d0
# ╠═80e407e4-f584-4253-b685-28e92b5c9119
# ╠═29da1d84-e6bc-4890-9a02-9c49e6505681
# ╠═9bee752e-9272-46af-8bc0-8c2a2a1807e2
# ╠═53a19766-1d66-43ee-b6df-3ef581746a78
# ╠═8c65263e-1e8e-40d8-b824-915b98a5ce57
# ╠═75f00c6e-2859-4b9e-8e6e-cb102f2fd071
# ╠═2486a3e8-6e80-4d84-a1b3-e8952f3f7faa
# ╠═571ede6a-ad6b-4710-a17d-9c3c7951d719
# ╠═a768ec77-2f9f-4e18-b488-b0cb756a4b00
# ╠═735efa70-2f53-45e0-a301-2ea716ea242a
# ╠═aea6c6c9-faae-4c7e-bf35-b1faace52177
# ╠═84389a19-1f13-4162-a51b-19dcf16ea0c1
# ╠═ee44f513-ffa0-4183-bcf9-fdbcb73dba69
# ╠═781bdb61-4311-4e5c-95d4-4deb751c7e6e
# ╠═aae30131-4654-4f98-83f1-9d0818189785
# ╠═6781ca9d-85fc-46c9-8059-e62e13edabd2
# ╠═ed515a86-9939-4123-9237-7747817b03d0
# ╠═b222c039-9c93-4bcb-a4cb-bc4219b11f8f
# ╠═392b2a08-597f-47e3-aeba-b45badcd43f5
# ╠═a44fce13-f88e-456b-b870-983503331d16
# ╠═341b5e17-056c-478d-8c3e-264bbf4bfeb0
# ╠═1294fa8e-d76b-406c-8d07-21e3aac0c59c
# ╠═6e36ec71-a706-4dab-be55-f2082f4a21c3
# ╠═0afc5cf7-852d-4d05-b3ad-7477bf86c7b9
# ╠═7f93a86c-17e4-4472-a9db-c6b43fce0dca
# ╠═4ff40840-74cb-4ac4-80b7-982c57478969
# ╠═600a0427-ad61-4317-aed9-41ba609c40c0
# ╠═57444eac-7b7b-47c5-8e85-1ca01b4230ce
# ╠═32abd57f-71b4-4188-9f28-0497c5ab7764
# ╠═6508d862-75c5-4076-9aec-29043a77493e
# ╠═daa2c7e2-2e6a-4a37-b32d-8569f21800bf
# ╠═5f1bb0e2-06e9-4476-8754-6175ed708b2e
# ╠═a7484a8b-6feb-4ccc-987b-d87f800970f2
# ╠═430ded66-dc00-411a-9db7-130f2d6555bf
# ╠═e1b3f8ef-9b44-4bec-8610-01145774be1a
# ╠═b0b44ab1-bd99-4239-b93f-3dea79d0c053
# ╠═a94baef6-c552-4b33-9139-0dc61dcf5aa9
# ╠═3cfc5665-936a-4621-8aa4-820c9a3e2c78
# ╠═2e606abc-4ae4-4798-88ff-48153e41038b
# ╠═d1b68ed8-f33e-4a02-b84f-03130ab96e58
# ╠═52a1f8c4-013a-4215-9b12-98f2dc23b2f3
# ╠═f2cd1fd6-3b73-4615-9d3f-078dcdf50911
# ╠═d51d3406-0e89-4a36-9591-199e3675d10f
# ╠═b9badf2d-b58e-4b12-919c-feb79f911926
# ╠═d192c0af-8551-4cba-b3c3-e8ea98c1b2a3
# ╠═cb4d2651-8471-4fe6-ae96-413306f9ecc8
# ╠═1639c8dc-e96d-4225-bcda-b838ad391d8b
# ╠═67ee804e-8877-4fc1-b2b5-a213bce545cd
# ╠═7e76da0f-91ca-4c5b-bf1d-ca10a2a6c5a8
# ╠═bfe80af1-a0ad-4702-9571-869aec104c28
# ╠═b5b6abab-40b0-430e-a012-59af58325a4d
# ╠═400df1dd-979a-45f4-9884-a6d9c012bd32
# ╠═f6bb72e2-267f-45e6-b2f9-905bda307b06
# ╠═dd49243c-3bbb-40fa-ac72-308c4ddebb57
# ╠═afa683f5-4bc3-4bfd-8e73-ba3e4c911368

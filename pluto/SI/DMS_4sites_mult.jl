### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 96c3c693-64d4-4380-9542-282b579a786a
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 443c5b5e-d5f2-48b4-8742-b7912595fcc3
using BioSequences: LongRNA

# ╔═╡ d87de269-6fe6-4612-bd31-a1413691fdf6
using Distributions: Gamma

# ╔═╡ 1acd1199-bc69-4038-b0ed-671d81233d32
using Makie: @L_str

# ╔═╡ ab006e15-451b-449f-b13e-ff011c647602
using NaNStatistics: nanmean, nanstd, nansum

# ╔═╡ fa4b7221-c125-4859-b8c9-d3a634154cf9
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 3a087371-61f0-4f30-8a01-895cb2a02f44
using Statistics: cor, mean

# ╔═╡ b1efaa23-f27c-487c-b815-6a6881f40116
md"# Imports"

# ╔═╡ 82f54c48-96d0-4cb4-b165-acfda769a05d
import PlutoUI, Makie, CairoMakie, StatsBase, SamApp2025

# ╔═╡ e9784767-7302-4860-ac5a-869fdceb0b14
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ a96298f7-d8ee-4ab7-9130-6c53fa8e4f53
PlutoUI.TableOfContents()

# ╔═╡ d84c3df6-991c-4f37-b748-26436f02ce92
md"# Load data"

# ╔═╡ 7011ccb1-7180-4a8f-a617-4ce5b60455cc
# load SHAPE data from Rep. 0, to get names of aptamers in DMS
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ d2ddfa18-5f7b-4308-a64b-da49815916de
# load SHAPE data for 500 aptamers, to get names of aptamers in DMS
shape_data_500 = SamApp2025.load_shapemapper_data_500v2_20240315();

# ╔═╡ 3997fa3e-3510-456a-aeb1-7e022f37f213
dms_df = SamApp2025.load_dms_data_sequences_table_20250303_with_aligned_sequences();

# ╔═╡ 0f697487-a276-4d63-b8fb-c386883fb822
dms_data = SamApp2025.load_dms_data_20250303();

# ╔═╡ be06284d-8899-4292-952e-985601efbd56
aptamer_names_500 = "APSAM-S2-" .* lpad.(0:499, 3, '0');

# ╔═╡ b8d90b73-f0f6-41f7-9a93-1d5e4b2002b3
dms_aptamer_origin = [
	begin
		if n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_rbm"]
			"rbm"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_seed70"]
			"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_full30"]
				"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_inf"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "rbm"]
			"rbm"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "infernal"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "Infrared"]
			"Infrared"
		end
	end for n = dms_data.aptamer_names
]

# ╔═╡ 39e55f6c-fc13-445c-9395-ee5edc521afb
dms_num_sites, dms_num_seqs, dms_num_conds = size(dms_data.shape_reactivities);

# ╔═╡ 9d43fb39-cf73-4a09-a38f-9577bb2db794
cond_sam = 1; cond_mg = 2;

# ╔═╡ f091912c-946b-4b0f-a3f4-b059c4a2fef8
nat_seqs = findall(dms_aptamer_origin .== "natural")

# ╔═╡ de8fe006-1990-4a3f-a52f-99a5d8fb5dbd
(; bps, nps, pks) = SamApp2025.RF00162_sites_paired()

# ╔═╡ 3474a86a-45e1-4af3-938e-9a463fc52887
all_dms_reactivities = [dms_data.shape_reactivities[i,n,cond_sam] for i=1:dms_num_sites for n=nat_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ 177f6584-1ec7-41c5-9f80-35195b523c61
bps_dms_reactivities = [dms_data.shape_reactivities[i,n,cond_sam] for i=bps for n=nat_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ e39daab1-e34f-43b1-ac41-7722f6fd7df5
nps_dms_reactivities = [dms_data.shape_reactivities[i,n,cond_sam] for i=nps for n=nat_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ f5a556fd-6da0-42f9-8a97-d6d9a8c00139
dms_stats = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = dms_data,
    paired_reactivities = bps_dms_reactivities,
    unpaired_reactivities = nps_dms_reactivities,
    all_reactivities = all_dms_reactivities,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 2c8f10b2-a588-49c3-a7d7-cdb2bd3ad864
md"# Paired and unpaired histogram datasets"

# ╔═╡ 0d65ddf8-6ea1-4c3d-87e7-21b088e7b70c
dms_log_odds_bp_sam = filter(isfinite, dms_stats.shape_log_odds[bps, nat_seqs, cond_sam])

# ╔═╡ cd51fd9f-722b-4eea-af4c-894b44ff6232
dms_log_odds_np_sam = filter(isfinite, dms_stats.shape_log_odds[nps, nat_seqs, cond_sam])

# ╔═╡ 2ea16128-b568-44dd-bfdb-0e0b5e3b2796
bps_reactivities_flat = filter(isfinite, dms_stats.shape_reactivities[bps, nat_seqs, cond_sam])

# ╔═╡ 064fef42-97ce-4de5-8c33-ec0227bcb47f
nps_reactivities_flat = filter(isfinite, dms_stats.shape_reactivities[nps, nat_seqs, cond_sam])

# ╔═╡ 4884e661-79ba-4720-aece-502304621b27
md"# Plots"

# ╔═╡ 6a315bbc-30a4-4760-a109-0ecd12227565
let fig = Makie.Figure()
	_sz = 150
	for (col, num_sites) = enumerate([1, 8, 24])

		data_bp = dropdims(sum(rand(dms_log_odds_bp_sam, num_sites, 10000); dims=1); dims=1)
		data_np = dropdims(sum(rand(dms_log_odds_np_sam, num_sites, 10000); dims=1); dims=1)

		bins = sqrt(num_sites) * (-20:0.05:10)

		ax = Makie.Axis(
			fig[1, col], width=_sz, height=_sz,
			xlabel = num_sites == 1 ? "Protection score" : "Total protection score",
			ylabel="Frequency", xgridvisible=false, ygridvisible=false,
			title= num_sites == 1 ? "1 site" : "$num_sites sites"
		)
		Makie.hist!(ax, data_bp; normalization=:pdf, bins, color=(:blue, 0.3))
		Makie.hist!(ax, data_np; normalization=:pdf, bins, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp; normalization=:pdf, bins, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np; normalization=:pdf, bins, color=:red, linewidth=2, label="unpaired")

		if num_sites == 1
			Makie.xlims!(ax, -3, 2)
			Makie.ylims!(ax, 0, 1)
			Makie.axislegend(ax; position=:rt, framevisible=false)
		elseif num_sites == 8
			Makie.xlims!(ax, -15, 10)
			Makie.ylims!(ax, 0, 0.2)
		else
			Makie.xlims!(ax, -30, 20)
			Makie.ylims!(ax, 0, 0.12)
		end

		####
		# Second row
		####

		data_bp = dropdims(sum(rand(bps_reactivities_flat, num_sites, 10000); dims=1); dims=1)
		data_np = dropdims(sum(rand(nps_reactivities_flat, num_sites, 10000); dims=1); dims=1)

		bins = sqrt(num_sites) * (-5:0.05:20)

		ax = Makie.Axis(
			fig[2, col], width=_sz, height=_sz,
			xlabel = num_sites == 1 ? "Reactivity" : "Total Reactivity",
			ylabel="Frequency", xgridvisible=false, ygridvisible=false,
			title= num_sites == 1 ? "1 site" : "$num_sites sites"
		)

		Makie.hist!(ax, data_bp; normalization=:pdf, bins, color=(:blue, 0.3))
		Makie.hist!(ax, data_np; normalization=:pdf, bins, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp; normalization=:pdf, bins, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np; normalization=:pdf, bins, color=:red, linewidth=2, label="unpaired")
		# Makie.stephist!(ax,
		# 	[dropdims(sum(rand(shape_log_odds_bp_sam, 8, 10000); dims=1); dims=1); dropdims(sum(rand(shape_log_odds_np_sam, 8, 10000); dims=1); dims=1)],
		# 	normalization=:pdf, bins=-11:0.05:5, color=:black, gap=-0.01
		# )
		if num_sites == 1
			Makie.xlims!(ax, -1, 2)
			Makie.ylims!(ax, 0, 5.8)
		elseif num_sites == 8
			Makie.xlims!(ax, -3, 10)
			Makie.ylims!(ax, 0, 0.6)
		else
			Makie.xlims!(ax, -3, 15)
			Makie.ylims!(ax, 0, 0.3)
		end
	end

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,2][1, 1, Makie.TopLeft()], "B)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,3][1, 1, Makie.TopLeft()], "C)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "D)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,2][1, 1, Makie.TopLeft()], "E)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,3][1, 1, Makie.TopLeft()], "F)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	save_path = mktempdir() * "/fig.pdf"
	Makie.save(save_path, fig)
	println("Saved figure to $save_path")
	fig
end

# ╔═╡ Cell order:
# ╠═b1efaa23-f27c-487c-b815-6a6881f40116
# ╠═96c3c693-64d4-4380-9542-282b579a786a
# ╠═82f54c48-96d0-4cb4-b165-acfda769a05d
# ╠═e9784767-7302-4860-ac5a-869fdceb0b14
# ╠═443c5b5e-d5f2-48b4-8742-b7912595fcc3
# ╠═d87de269-6fe6-4612-bd31-a1413691fdf6
# ╠═1acd1199-bc69-4038-b0ed-671d81233d32
# ╠═ab006e15-451b-449f-b13e-ff011c647602
# ╠═fa4b7221-c125-4859-b8c9-d3a634154cf9
# ╠═3a087371-61f0-4f30-8a01-895cb2a02f44
# ╠═a96298f7-d8ee-4ab7-9130-6c53fa8e4f53
# ╠═d84c3df6-991c-4f37-b748-26436f02ce92
# ╠═7011ccb1-7180-4a8f-a617-4ce5b60455cc
# ╠═d2ddfa18-5f7b-4308-a64b-da49815916de
# ╠═3997fa3e-3510-456a-aeb1-7e022f37f213
# ╠═0f697487-a276-4d63-b8fb-c386883fb822
# ╠═be06284d-8899-4292-952e-985601efbd56
# ╠═b8d90b73-f0f6-41f7-9a93-1d5e4b2002b3
# ╠═39e55f6c-fc13-445c-9395-ee5edc521afb
# ╠═9d43fb39-cf73-4a09-a38f-9577bb2db794
# ╠═f091912c-946b-4b0f-a3f4-b059c4a2fef8
# ╠═de8fe006-1990-4a3f-a52f-99a5d8fb5dbd
# ╠═3474a86a-45e1-4af3-938e-9a463fc52887
# ╠═177f6584-1ec7-41c5-9f80-35195b523c61
# ╠═e39daab1-e34f-43b1-ac41-7722f6fd7df5
# ╠═f5a556fd-6da0-42f9-8a97-d6d9a8c00139
# ╠═2c8f10b2-a588-49c3-a7d7-cdb2bd3ad864
# ╠═0d65ddf8-6ea1-4c3d-87e7-21b088e7b70c
# ╠═cd51fd9f-722b-4eea-af4c-894b44ff6232
# ╠═2ea16128-b568-44dd-bfdb-0e0b5e3b2796
# ╠═064fef42-97ce-4de5-8c33-ec0227bcb47f
# ╠═4884e661-79ba-4720-aece-502304621b27
# ╠═6a315bbc-30a4-4760-a109-0ecd12227565

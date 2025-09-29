### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 1a0cd11b-470f-4cdc-9717-c06bb1b809a0
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ b12e4a16-d0e0-4348-aaa0-8d8f10c0a58e
using BioSequences: LongRNA

# ╔═╡ 5af001ab-2b23-44ac-9460-5fffc438ff4e
using Distributions: Gamma

# ╔═╡ 52b29ec7-4980-4edc-a430-384c083e6352
using Makie: @L_str

# ╔═╡ b4702133-da4c-405d-b163-78f6abb64d02
using NaNStatistics: nanmean

# ╔═╡ 443a0c3e-19a6-46d0-92ed-9bcc6a51b9a6
using NaNStatistics: nanstd

# ╔═╡ b1091fdf-b182-4cd0-b4bb-9c54ac910772
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 21f90a50-a26f-4a95-94cd-fcdefca4e4ff
using Statistics: cor

# ╔═╡ e3210ca1-214f-4c20-9daa-ae0b1fc359d0
using Statistics: mean

# ╔═╡ dade2a10-73a2-11ef-1f86-170e1e6ad9d2
md"# Imports"

# ╔═╡ c566a92b-00a9-4e12-913c-9969ea085d46
import PlutoUI

# ╔═╡ 90273f59-a2d8-4d88-a131-2e2e713dd17a
import CairoMakie

# ╔═╡ 92330187-3f1e-491c-a322-e9ea725094b4
import Makie

# ╔═╡ 37626325-b59b-492d-8808-5f8f5f5b6f6f
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 389cfe17-9459-47eb-8163-fa1f3e050814
import SamApp2025

# ╔═╡ 39dc0020-0f5c-49de-8c2e-4e0da6f7e72e
import StatsBase

# ╔═╡ ada2e001-3480-43e5-82a1-d433257446ce
PlutoUI.TableOfContents()

# ╔═╡ 0220c511-1c83-464c-8c0e-23225ae05b11
md"# Load data"

# ╔═╡ 1f0e0eca-913a-481c-b0a5-a23f3bf25609
# load SHAPE data
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ ba387702-7fdb-4e1b-97d3-644c22962191
# split rep0 from rep4+5
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 58ae7a5f-6660-44f2-bbc2-bea8dfc2be96
# split rep0 from rep4+5
shape_data_rep45 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 54c7db92-3427-48b7-99c0-eb9f95e283de
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ c8cd8fda-5472-4a86-8ac7-951c4028dd3a
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ e51867da-1047-4d88-976b-4073244b8e68
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 75219269-12be-48cf-b757-c8859b1b60b3
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 34011e0b-c892-479b-94d6-644976ac7b44
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ e8a35f14-f687-41d1-8dd0-414d7ff0a8f2
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 8839245b-3973-4220-8c29-7e3b92a6312d
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ ad0abe4d-a689-4ba0-a5b9-a2b782ef39d7
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ 0bf4f3b4-d3ab-4dab-ba2c-e85123007ce1
(; bps, nps, pks) = SamApp2025.RF00162_sites_paired()

# ╔═╡ 5365c629-5c78-4815-897a-f074c39b9f89
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ a6ba0af0-0c96-4344-8b75-c10bc003dc88
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 235ace16-e96a-48b7-90cb-ecf101f57a30
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 643003a0-2a83-4f69-906a-2e31844f7d1c
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 0fc892a2-1d6a-4b4d-9442-7bfcdc077ed6
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ bb29a1f1-b3e3-4e7f-bc79-0d1555005483
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, seed_seqs, conds_sam_rep0];

# ╔═╡ 55f0f62e-696e-44c3-9098-75e3a0041f2b
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, seed_seqs, conds_sam_rep0];

# ╔═╡ 2851791f-c0ae-43c8-b69f-f8791a3fb756
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ 8e3703f3-d8b4-46fc-8b5c-be5a9a216d3c
shape_stats_rep0 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ d526be42-aca6-4dcb-a75a-428d1a318568
md"# Figure"

# ╔═╡ a3c5a2ca-f066-497b-86ce-764367485a97
# structural motifs
struct_bands = [
    (; x0=0.5, xf=8.5, color="blue", alpha=0.1), # P1
    (; x0=100.5, xf=108.5, color="blue", alpha=0.1), # P1
    (; x0=11.5, xf=16.5, color="green", alpha=0.1), # P2
    (; x0=20.5, xf=23.5, color="green", alpha=0.1), # P2
    (; x0=28.5, xf=31.5, color="green", alpha=0.1), # P2
    (; x0=37.5, xf=42.5, color="green", alpha=0.1), # P2
    (; x0=42.5, xf=46.5, color="orange", alpha=0.1), # P3
    (; x0=47.5, xf=53.5, color="orange", alpha=0.1), # P3
    (; x0=60.5, xf=64.5, color="orange", alpha=0.1), # P3
    (; x0=66.5, xf=72.5, color="orange", alpha=0.1), # P3
    (; x0=80.5, xf=86.5, color="teal", alpha=0.1), # P4
    (; x0=91.5, xf=97.5, color="teal", alpha=0.1), # P4
    (; x0=24.5, xf=28.5, color="red", alpha=0.1), # Pk
    (; x0=76.5, xf=80.5, color="red", alpha=0.1), # Pk
];

# ╔═╡ 22ef8b87-23f2-415b-8a8a-7ae8cd001468
function __ss_bands!(ax)
    for (x0, xf, color, alpha) = struct_bands
        #Makie.band!(ax, x0:xf, 0, 1, color=(:lightblue, 0.5))
        Makie.vspan!(ax, x0, xf; color=(color, alpha), ymin=0.0, ymax=1)
    end

    # helices = [
    #     1:8, 102:108, # P1
    #     13:17, 21:23, 29:31, 38:42, # P2
    #     43:46, 48:53, 61:64, 67:72, # P3
    #     81:86, 92:98 # P4
    # ]
    # loops = [
    #     9:12, 73:76, # Central loop
    #     18:20, 32:37, # L2
    #     54:60, # L3
    #     87:91, 99:101, # L4
    # ]

    # pk = [24:28, 77:80]
    # for ps in helices
    #     Makie.band!(ax, (first(ps) - 0.5):(last(ps) + 0.5), 0, 1, color=(:lightblue, 0.5))
    # end
    # for ps in loops
    #     Makie.band!(ax, (first(ps) - 0.5):(last(ps) + 0.5), -2, 0, color=(:red, 0.3))
    # end
    # for ps in pk
    #     Makie.band!(ax, (first(ps) - 0.5):(last(ps) + 0.5), -2, 1, color=(:orange, 0.3))
    # end
    # bulges
    # Makie.band!(ax, (47 - 0.5):(47 + 0.5), -2, 1, color=(:cyan, 0.3))
    # Makie.band!(ax, (65 - 0.5):(66 + 0.5), -2, 1, color=(:cyan, 0.3))
end

# ╔═╡ 6d20a823-1c34-4d10-9404-8e9254470682
_thresh = 0.1

# ╔═╡ fe297c97-3184-4ff7-8ca4-26ecb7e0b590
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1], width=800, height=120, xlabel="site", ylabel="protection score", xticks=5:10:108, yticks=[-1.5, -0.5, 0.5], xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	#Makie.hlines!(ax, 0, color=(:silver, 0.5))
	for (x0, xf, color, alpha) = struct_bands
	    #Makie.poly!(ax, Makie.Point2f[(x0, 0), (xf, 0), (xf, 1), (x0, 1)]; color=(color, alpha))
	    Makie.vspan!(ax, x0, xf; color=(color, alpha))
	end

	S_mg = nanmean(shape_stats_rep0.shape_log_odds[:, nat_seqs, conds_mg_rep0]; dim=(2,3))
	S_sam = nanmean(shape_stats_rep0.shape_log_odds[:, nat_seqs, conds_sam_rep0]; dim=(2,3))

	S_mg_err = nanstd(shape_stats_rep0.shape_log_odds[:, nat_seqs, conds_mg_rep0]; dim=(2,3))
	S_sam_err = nanstd(shape_stats_rep0.shape_log_odds[:, nat_seqs, conds_sam_rep0]; dim=(2,3))

	Makie.hlines!(ax, 0.0; linestyle=:dash, color=(:red, 0.5))
	Makie.stairs!(ax, 1:108, S_mg; step=:center, color=:gray)
	Makie.stairs!(ax, 1:108, S_sam; step=:center, color=:purple)
	Makie.errorbars!(ax, 1:108, S_mg, S_mg_err/2, color=(:gray, 0.5), linewidth=1)
	Makie.errorbars!(ax, 1:108, S_sam, S_sam_err/2, color=(:purple, 0.5), linewidth=1)
	Makie.xlims!(ax, 0, 109)
	Makie.ylims!(ax, -1.5, 0.7)
	Makie.hidespines!(ax, :t, :r)

	# _thresh = 0.1
	# #_sites = findall((R_mg .< -_thresh) .& (R_sam .> _thresh))
	_sites = SamApp2025.hallmark_sites_20230507
	Makie.scatter!(ax, _sites, -1.4one.(_sites), markersize=10, color=:black, marker=:utriangle)


	ax = Makie.Axis(fig[2,1], width=800, height=120, xlabel="site", ylabel="Reactivity", xticks=5:10:108, yticks=-1:5, xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    #Makie.poly!(ax, Makie.Point2f[(x0, 0), (xf, 0), (xf, 1), (x0, 1)]; color=(color, alpha))
	    Makie.vspan!(ax, x0, xf; color=(color, alpha))
	end

	R_mg = nanmean(shape_stats_rep0.shape_reactivities[:, nat_seqs, conds_mg_rep0]; dim=(2,3))
	R_sam = nanmean(shape_stats_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0]; dim=(2,3))

	R_mg_err = nanstd(shape_stats_rep0.shape_reactivities[:, nat_seqs, conds_mg_rep0]; dim=(2,3))
	R_sam_err = nanstd(shape_stats_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0]; dim=(2,3))

	Makie.stairs!(ax, 1:108, R_mg; step=:center, color=:gray)
	Makie.stairs!(ax, 1:108, R_sam; step=:center, color=:purple)
	Makie.errorbars!(ax, 1:108, R_mg, R_mg_err/2, color=(:gray, 0.5), linewidth=1)
	Makie.errorbars!(ax, 1:108, R_sam, R_sam_err/2, color=(:purple, 0.5), linewidth=1)
	Makie.scatter!(ax, _sites, -0.5one.(_sites), markersize=10, color=:black, marker=:utriangle)
	Makie.xlims!(ax, 0  , 109)
	#Makie.ylims!(ax, -1.5, 0.7)
	Makie.hidespines!(ax, :t, :r)

	# _thresh = 0.1
	# #_sites = findall((R_mg .< -_thresh) .& (R_sam .> _thresh))
	_sites = SamApp2025.hallmark_sites_20230507
	#Makie.scatter!(ax, _sites, -1.4one.(_sites), markersize=10, color=:green, marker=:utriangle)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Hallmark_sites.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═dade2a10-73a2-11ef-1f86-170e1e6ad9d2
# ╠═1a0cd11b-470f-4cdc-9717-c06bb1b809a0
# ╠═c566a92b-00a9-4e12-913c-9969ea085d46
# ╠═90273f59-a2d8-4d88-a131-2e2e713dd17a
# ╠═92330187-3f1e-491c-a322-e9ea725094b4
# ╠═37626325-b59b-492d-8808-5f8f5f5b6f6f
# ╠═389cfe17-9459-47eb-8163-fa1f3e050814
# ╠═39dc0020-0f5c-49de-8c2e-4e0da6f7e72e
# ╠═b12e4a16-d0e0-4348-aaa0-8d8f10c0a58e
# ╠═5af001ab-2b23-44ac-9460-5fffc438ff4e
# ╠═52b29ec7-4980-4edc-a430-384c083e6352
# ╠═b4702133-da4c-405d-b163-78f6abb64d02
# ╠═443a0c3e-19a6-46d0-92ed-9bcc6a51b9a6
# ╠═b1091fdf-b182-4cd0-b4bb-9c54ac910772
# ╠═21f90a50-a26f-4a95-94cd-fcdefca4e4ff
# ╠═e3210ca1-214f-4c20-9daa-ae0b1fc359d0
# ╠═ada2e001-3480-43e5-82a1-d433257446ce
# ╠═0220c511-1c83-464c-8c0e-23225ae05b11
# ╠═1f0e0eca-913a-481c-b0a5-a23f3bf25609
# ╠═ba387702-7fdb-4e1b-97d3-644c22962191
# ╠═58ae7a5f-6660-44f2-bbc2-bea8dfc2be96
# ╠═54c7db92-3427-48b7-99c0-eb9f95e283de
# ╠═c8cd8fda-5472-4a86-8ac7-951c4028dd3a
# ╠═e51867da-1047-4d88-976b-4073244b8e68
# ╠═75219269-12be-48cf-b757-c8859b1b60b3
# ╠═34011e0b-c892-479b-94d6-644976ac7b44
# ╠═e8a35f14-f687-41d1-8dd0-414d7ff0a8f2
# ╠═8839245b-3973-4220-8c29-7e3b92a6312d
# ╠═ad0abe4d-a689-4ba0-a5b9-a2b782ef39d7
# ╠═0bf4f3b4-d3ab-4dab-ba2c-e85123007ce1
# ╠═5365c629-5c78-4815-897a-f074c39b9f89
# ╠═a6ba0af0-0c96-4344-8b75-c10bc003dc88
# ╠═235ace16-e96a-48b7-90cb-ecf101f57a30
# ╠═643003a0-2a83-4f69-906a-2e31844f7d1c
# ╠═0fc892a2-1d6a-4b4d-9442-7bfcdc077ed6
# ╠═bb29a1f1-b3e3-4e7f-bc79-0d1555005483
# ╠═55f0f62e-696e-44c3-9098-75e3a0041f2b
# ╠═2851791f-c0ae-43c8-b69f-f8791a3fb756
# ╠═8e3703f3-d8b4-46fc-8b5c-be5a9a216d3c
# ╠═d526be42-aca6-4dcb-a75a-428d1a318568
# ╠═a3c5a2ca-f066-497b-86ce-764367485a97
# ╠═22ef8b87-23f2-415b-8a8a-7ae8cd001468
# ╠═6d20a823-1c34-4d10-9404-8e9254470682
# ╠═fe297c97-3184-4ff7-8ca4-26ecb7e0b590

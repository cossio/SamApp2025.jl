### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ b3c73892-1482-11ef-1248-f70b8e38afbd
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 38405504-270a-4b53-83f9-c250f84430e7
using BioSequences: LongRNA

# ╔═╡ 53290a43-f19a-43e0-a4be-7be73dcc339b
using BioSequences: RNA_Gap

# ╔═╡ 4bc1fbc6-c7b3-471d-bd86-6f7a59cc0749
using SamApp2025: onehot

# ╔═╡ f2391009-c7cf-4ea2-a7ca-cd9f72605d32
using DataFrames: DataFrame

# ╔═╡ 4d57ef03-6516-4e43-908b-d31f045ec59b
using LinearAlgebra: eigen

# ╔═╡ 7f08a116-eae3-4fd0-b395-d13038997465
using NaNStatistics: nanmean, nanstd, nansum

# ╔═╡ e162b72f-d324-4725-abc4-c8ab73770cf1
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ e908abe1-dd9e-4f5d-bb81-941d6edd08d9
using Statistics: cor

# ╔═╡ d806c591-00b9-49ba-a755-e3923ae41670
using Statistics: mean

# ╔═╡ 0a5f759b-c834-4a5b-899c-0988c03cb757
using SamApp2025: successes_tuple_str

# ╔═╡ ca2287c5-2509-4870-a335-605c7ec1022d
md"# Imports"

# ╔═╡ 60221493-44cb-41ff-b871-3352d42b0c21
import CairoMakie

# ╔═╡ 415f361f-8ae0-40d7-93d6-2f0f1ef2ee4a
import CSV

# ╔═╡ 5a94fc1c-d82c-4bb0-b8be-3c28c13e0254
import Dates

# ╔═╡ 517a5648-ffb4-401d-b561-60a7fe49db11
import FASTX

# ╔═╡ ae442ee6-cccc-4447-8763-14462257a882
import Infernal

# ╔═╡ 845a563a-1853-4c00-b5ed-7768d8f852b5
import Makie

# ╔═╡ a7379a53-a3e5-4e8e-9c9f-1fa9c58ff8a1
import Rfam

# ╔═╡ ac311ee6-062d-4e88-9c58-8852126c0924
import SamApp2025

# ╔═╡ d354f10d-2eb1-4e74-bcb9-a4f01028cdf7
import PlutoUI

# ╔═╡ b7e64deb-3140-4038-abaf-6a7d2996358c
PlutoUI.TableOfContents()

# ╔═╡ 85110743-90c3-40a1-ab08-abc481480223
md"""
# Load data
"""

# ╔═╡ b214a89f-7a82-46ee-8b04-2a55cc257bc9
shape_data = SamApp2025.load_shapemapper_data_500v2_20240315();

# ╔═╡ 07bdfaa2-587f-49b3-a50d-e92618399eb1
RF00162_hits = SamApp2025.rfam_RF00162_hits();

# ╔═╡ 41636ccb-731a-453e-b4cc-f5643de8a63d
rbm_energies = free_energy(SamApp2025.rbm2022(), onehot(shape_data.aligned_sequences));

# ╔═╡ 8b8d735b-dd49-4727-915b-795ef7a987ae
aptamer_natural_distances = SamApp2025.hamming(onehot(shape_data.aligned_sequences), onehot(RF00162_hits));

# ╔═╡ 8f39ec27-bdff-4762-81c1-0d869074ff54
conds_sam = [1,2];

# ╔═╡ c8187a90-68aa-4e7f-af58-5765f37d88e1
conds_mg = [4];

# ╔═╡ 161737a5-1f5b-4497-b00f-afc85da70569
conds_30C = [6];

# ╔═╡ 7252e0bc-c0cb-4477-a108-0896ae66c7ea
the_conds = vcat(conds_sam, conds_mg, conds_30C)

# ╔═╡ 6220c7db-be0c-47b6-b59e-748eabaadd87
shape_data.conditions[the_conds]

# ╔═╡ adaec625-d75e-429d-b6a2-18ec51c80c4a
bps, nps, pks = SamApp2025.RF00162_sites_paired()

# ╔═╡ dd4b8dd3-f151-4564-af94-8102cbf6917f
ss_sites = SamApp2025.RF00162_sites_annotated_secondary_structure()

# ╔═╡ 5604af2b-1435-4c64-92d8-546e758814bd
bps_reactivities = shape_data.shape_reactivities[bps, :, conds_sam];

# ╔═╡ 99b96798-4a5a-4f75-a576-b03488487e02
nps_reactivities = shape_data.shape_reactivities[nps, :, conds_sam];

# ╔═╡ b756b05f-b4f5-4a64-b344-bf360e9a8bd3
all_reactivities = shape_data.shape_reactivities[:, :, conds_sam];

# ╔═╡ abe07005-c0b9-45d2-869c-680af950bf5f
shape_stats = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data,
    paired_reactivities = bps_reactivities,
    unpaired_reactivities = nps_reactivities,
    all_reactivities = all_reactivities,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ e4dc6b05-8f9d-4177-87a8-7c8ae622894b
_sites = SamApp2025.hallmark_sites_20230507

# ╔═╡ f7517c96-ae0d-43bb-9f12-f0d077172e70
p4_pos = union(81:86, 92:97); # P4 helix positions

# ╔═╡ 76bfd0d4-f723-44d8-b2a5-3f2645a9e8cb
# length of P4 segment for probed sequences
aptamer_p4_length = length(p4_pos) .- [count(seq[p4_pos] .== RNA_Gap) for seq = shape_data.aligned_sequences];

# ╔═╡ 748afc4e-a8f7-4a12-99aa-296aff325520
x_mg = nansum(shape_stats.shape_log_odds[_sites, :,  conds_mg]; dim=(1,3));

# ╔═╡ 96b195c0-685e-48ec-a9f2-92e69e17b31c
x_sam = nansum(shape_stats.shape_log_odds[_sites, :, conds_sam]; dim=(1,3));

# ╔═╡ 4f3ed2b1-4089-4caa-8180-d584b739d66f
_thresh = log(5)

# ╔═╡ 8c528821-ba3f-45e2-965d-e93aab90fc72
_responds_sam_yes = (x_mg .< -_thresh) .& (x_sam .> +_thresh);

# ╔═╡ b33cb49f-df54-44e9-af9e-bd2b6aa1e42a
_responds_sam_nop = (x_mg .> +_thresh) .| (x_sam .< -_thresh);

# ╔═╡ 762b903f-0ded-4e3b-804f-59fc73c37dfa
df = DataFrame(;
    aligned_sequences = shape_data.aligned_sequences,
    aptamer_origin = shape_data.aptamer_origin,
    switch_yes = _responds_sam_yes,
    switch_no = _responds_sam_nop
)

# ╔═╡ 96349783-915d-4328-bc41-4e8d7e516ad7
rbm_seqs = findall(shape_data.aptamer_origin .== "rbm");

# ╔═╡ 21966bd0-143a-4e95-9881-b645b9b5a261
inf_seqs = findall(shape_data.aptamer_origin .== "infernal");

# ╔═╡ 53653580-28c9-46ce-a9e8-63ec523e8022
inf_untangled_seqs = findall(shape_data.aptamer_criteria .== "perm");

# ╔═╡ b4a38a4a-45ae-4666-803b-0769c0ea06de
inf_uniform_seqs = findall(shape_data.aptamer_criteria .== "uniform");

# ╔═╡ 3ac42368-de6d-4012-98ab-cf275a9436a8
rbm_p4_seqs = findall(shape_data.aptamer_criteria .== "p4");

# ╔═╡ be2b3c00-d909-48e9-b009-59f2fd996534
# CM model from Rfam (this has the noisy floor!)
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162");

# ╔═╡ dd238bab-06c5-453a-9a12-bd0711b82c42
RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")

# ╔═╡ 81926f04-dd75-4e09-b1b2-e39efbf48151
RF00162_seed_match_cols = findall(≠('.'), SamApp2025.stockholm_ss(RF00162_seed_stk.out));

# ╔═╡ 111eca50-bfdd-42e7-b495-a4ee8cc6d48b
RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM") # WARNING: this has inserts marked as '-'

# ╔═╡ d63552db-2eed-46e0-ad6a-801644e6ca2a
RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))

# ╔═╡ ea18fef7-c1c6-46b9-b362-94e6762fc52f
RF00162_seed_seqs_noinserts = LongRNA{4}.([FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records]);

# ╔═╡ 54ca7e3e-8b56-4b77-9488-f0c4ea98a1e4
# trimmed (no inserts) aligned fasta
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");

# ╔═╡ 11478b83-6854-466c-8817-9fced3bc4b7b
# these are already aligned and without inserts
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ 1ee2688f-0fb4-4227-b101-156b310a34ee
# sites that have some non-zero fluctuations
# We need to separate frozen sites below because otherwise cor and eigen give NaN, infinities, and fail
_variable_sites_flag = vec(all(0 .< mean(SamApp2025.onehot(RF00162_hits_sequences); dims=3) .< 1; dims=1));

# ╔═╡ 49ca550c-5662-4ae3-841e-4411c670c6b7
_variable_sites = findall(_variable_sites_flag);

# ╔═╡ aab4cc7e-15c0-43a3-b9cb-35e1559478d1
RF00162_hits_var_sites_only = SamApp2025.onehot(RF00162_hits_sequences)[:, _variable_sites, :];

# ╔═╡ d82b6d02-a0ef-418b-898f-6dca245ff100
RF00162_hits_cor = cor(reshape(RF00162_hits_var_sites_only, :, size(RF00162_hits_var_sites_only, 3)); dims=2);

# ╔═╡ d2b92466-6365-41b7-822e-2bf2df895dd3
RF00162_hits_eig = eigen(RF00162_hits_cor);

# ╔═╡ dcf44d38-95df-4e27-b4ef-3858cb0a40ac
# remap the variable sites eigenvectors back to the original consensus sequence numbering
begin
	RF00162_hits_eigvec = zeros(5, 108, size(RF00162_hits_eig.vectors, 1))
	for n in 1:size(RF00162_hits_eig.vectors, 1)
	    vec(view(RF00162_hits_eigvec, :, _variable_sites, n)) .= RF00162_hits_eig.vectors[:, n]
	end
end

# ╔═╡ 2c4015ad-8778-44f9-9d27-58a5ccbb5bd5
shape_sequences_onehot = SamApp2025.onehot(LongRNA{4}.(shape_data.aligned_sequences));

# ╔═╡ b588e294-6449-45c7-afbc-6535f1dfa6eb
__proj_500 = reshape(shape_sequences_onehot, 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ e6028e8e-e58e-475e-994f-60bdecf164eb
__proj_hits = reshape(SamApp2025.onehot(RF00162_hits_sequences), 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ 91597c20-dace-4a3b-854a-8682173d9476
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

# ╔═╡ 57fea3b4-eb6c-4064-87b5-314d83b1e117
md"# Plots"

# ╔═╡ 36eb9427-db56-4265-9621-0a91efe99a93
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1][1,1], xlabel="PC1", ylabel="PC2", width=300, height=300, xgridvisible=false, ygridvisible=false, title="Probed sequences")
	Makie.scatter!(ax, __proj_hits[:, end], __proj_hits[:, end - 1], markersize=10, color=(:gray, 0.5), label="MSA")
	# Makie.scatter!(ax,
	#     __proj_500[(aptamer_p4_length .< 1) .& _responds_sam_yes, end],
	#     __proj_500[(aptamer_p4_length .< 1) .& _responds_sam_yes, end - 1],
	#     markersize=15, color=:orange, marker='●'
	# )
	Makie.scatter!(ax,
	    __proj_500[shape_data.aptamer_origin .== "rbm", end],
	    __proj_500[shape_data.aptamer_origin .== "rbm", end - 1],
	    markersize=10, color=:blue, label="RBM", marker=:cross
	)
	Makie.scatter!(ax,
	    __proj_500[inf_uniform_seqs, end],
	    __proj_500[inf_uniform_seqs, end - 1],
	    markersize=8, color=:red, label="dCM", marker=:cross
	)
	Makie.scatter!(ax,
	    __proj_500[inf_untangled_seqs, end],
	    __proj_500[inf_untangled_seqs, end - 1],
	    markersize=8, color=:orange, label="uCM", marker=:cross
	)
	Makie.axislegend(ax, position=:lb, framevisible=false, patchlabelgap=-3)

	ax = Makie.Axis(fig[1,1][2,1], width=300, height=300, xlabel="Divergence from closest natural", ylabel="RBM score", xticks=0:0.1:0.6, yticks=200:25:350, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances; dims=2)) / 108, -rbm_energies, markersize=15, color=(:teal, 0.2), marker='●', label="All probed")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[rbm_seqs ∩ findall(_responds_sam_yes), :]; dims=2)) / 108, -rbm_energies[rbm_seqs ∩ findall(_responds_sam_yes)], markersize=10, color=:blue, marker='●', label="RBM (✓)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[rbm_seqs ∩ findall(_responds_sam_nop), :]; dims=2)) / 108, -rbm_energies[rbm_seqs ∩ findall(_responds_sam_nop)], markersize=10, color=:blue, marker='O', label="RBM (❌)")
	# Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_seqs ∩ findall(_responds_sam_yes), :]; dims=2)), -rbm_energies[inf_seqs ∩ findall(_responds_sam_yes)], markersize=10, color=:red, marker='●', label="CM (✓)")
	# Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_seqs ∩ findall(_responds_sam_nop), :]; dims=2)), -rbm_energies[inf_seqs ∩ findall(_responds_sam_nop)], markersize=10, color=:red, marker='O', label="CM (❌)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_uniform_seqs ∩ findall(_responds_sam_yes), :]; dims=2)) / 108, -rbm_energies[inf_uniform_seqs ∩ findall(_responds_sam_yes)], markersize=10, color=:red, marker='●', label="dCM (✓)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_uniform_seqs ∩ findall(_responds_sam_nop), :]; dims=2)) / 108, -rbm_energies[inf_uniform_seqs ∩ findall(_responds_sam_nop)], markersize=10, color=:red, marker='O', label="dCM (❌)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_untangled_seqs ∩ findall(_responds_sam_yes), :]; dims=2)) / 108, -rbm_energies[inf_untangled_seqs ∩ findall(_responds_sam_yes)], markersize=10, color=:orange, marker='●', label="uCM (✓)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances[inf_untangled_seqs ∩ findall(_responds_sam_nop), :]; dims=2)) / 108, -rbm_energies[inf_untangled_seqs ∩ findall(_responds_sam_nop)], markersize=10, color=:orange, marker='O', label="uCM (❌)")
	Makie.ylims!(ax, 240, 360)
	Makie.axislegend(ax, position=(-0.02, -0.01), framevisible=false, nbanks=1, colgap=1, rowgap=0.1, patchlabelgap=0)


	_dummy_axis = Makie.Axis(fig[1,2][0,1])
	Makie.hidespines!(_dummy_axis, :t, :b, :l, :r)
	Makie.hidexdecorations!(_dummy_axis)
	Makie.hideydecorations!(_dummy_axis)

	_width = 500
	_height = 65

	# example of no P4
	n_ex = 116

	_R_mg = shape_data.shape_reactivities[:, n_ex, only(conds_mg)]
	_R_sam = shape_data.shape_reactivities[:, n_ex, conds_sam[2]]

	ax_react_1 = Makie.Axis(fig[1,2][1,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:4:8, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_1, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_1, 1:108, _R_mg, color=:gray, step=:center, label="no SAM")
	Makie.stairs!(ax_react_1, 1:108, _R_sam, color=:purple, step=:center, label="with SAM")
	Makie.hidespines!(ax_react_1, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_1)
	#Makie.axislegend(ax_react_1, position=(0.0, -13), framevisible=false)

	ax_diff_1 = Makie.Axis(fig[1,2][2,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-10:2:10, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_1, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_1, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff_1, _sites, -2.5one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff_1, 0, 109)
	Makie.hidespines!(ax_diff_1, :r, :b, :t)
	Makie.hidexdecorations!(ax_diff_1)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	# example distant
	n_ex = 284

	_R_mg = shape_data.shape_reactivities[:, n_ex, only(conds_mg)]
	_R_sam = shape_data.shape_reactivities[:, n_ex, conds_sam[2]]

	ax_react_2 = Makie.Axis(fig[1,2][3,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, yticks=0:4:8, ylabel="react.", xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_2, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_2, 1:108, _R_mg, color=:gray, step=:center, label="no SAM")
	Makie.stairs!(ax_react_2, 1:108, _R_sam, color=:purple, step=:center, label="with SAM")
	Makie.hidespines!(ax_react_2, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_2)

	ax_diff_2 = Makie.Axis(fig[1,2][4,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, yticks=-10:1:10, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_2, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_2, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff_2, _sites, -1.65one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.hidespines!(ax_diff_2, :r, :t)
	#Makie.scatter!(ax_diff_2, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react_1, ax_diff_1)
	Makie.linkxaxes!(ax_react_2, ax_diff_2)
	#Makie.linkyaxes!(ax_react_1, ax_react_2)
	#Makie.linkyaxes!(ax_diff_1, ax_diff_2)

	# Makie.ylims!(ax_diff_1, -1.5, 1)
	# Makie.ylims!(ax_diff_2, -1.5, 1)
	Makie.ylims!(ax_diff_1, -2.7, 3.4)
	Makie.ylims!(ax_diff_2, -1.8, 0.7)

	Makie.ylims!(ax_react_1, -0.5, 8)
	Makie.ylims!(ax_react_2, -0.5, 8)

	Makie.xlims!(ax_react_1, 0.5, 108.5)
	Makie.xlims!(ax_react_2, 0.5, 108.5)
	Makie.xlims!(ax_diff_1,  0.5, 108.5)
	Makie.xlims!(ax_diff_2,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/500 aptamers.pdf", fig)
	save_path = mktempdir() * "/fig.pdf"
	println("Saving to $save_path")
	Makie.save(save_path, fig)
	fig
end

# ╔═╡ Cell order:
# ╠═ca2287c5-2509-4870-a335-605c7ec1022d
# ╠═b3c73892-1482-11ef-1248-f70b8e38afbd
# ╠═60221493-44cb-41ff-b871-3352d42b0c21
# ╠═415f361f-8ae0-40d7-93d6-2f0f1ef2ee4a
# ╠═5a94fc1c-d82c-4bb0-b8be-3c28c13e0254
# ╠═517a5648-ffb4-401d-b561-60a7fe49db11
# ╠═ae442ee6-cccc-4447-8763-14462257a882
# ╠═845a563a-1853-4c00-b5ed-7768d8f852b5
# ╠═a7379a53-a3e5-4e8e-9c9f-1fa9c58ff8a1
# ╠═ac311ee6-062d-4e88-9c58-8852126c0924
# ╠═d354f10d-2eb1-4e74-bcb9-a4f01028cdf7
# ╠═38405504-270a-4b53-83f9-c250f84430e7
# ╠═53290a43-f19a-43e0-a4be-7be73dcc339b
# ╠═4bc1fbc6-c7b3-471d-bd86-6f7a59cc0749
# ╠═f2391009-c7cf-4ea2-a7ca-cd9f72605d32
# ╠═4d57ef03-6516-4e43-908b-d31f045ec59b
# ╠═7f08a116-eae3-4fd0-b395-d13038997465
# ╠═e162b72f-d324-4725-abc4-c8ab73770cf1
# ╠═e908abe1-dd9e-4f5d-bb81-941d6edd08d9
# ╠═d806c591-00b9-49ba-a755-e3923ae41670
# ╠═0a5f759b-c834-4a5b-899c-0988c03cb757
# ╠═b7e64deb-3140-4038-abaf-6a7d2996358c
# ╠═85110743-90c3-40a1-ab08-abc481480223
# ╠═b214a89f-7a82-46ee-8b04-2a55cc257bc9
# ╠═07bdfaa2-587f-49b3-a50d-e92618399eb1
# ╠═41636ccb-731a-453e-b4cc-f5643de8a63d
# ╠═8b8d735b-dd49-4727-915b-795ef7a987ae
# ╠═8f39ec27-bdff-4762-81c1-0d869074ff54
# ╠═c8187a90-68aa-4e7f-af58-5765f37d88e1
# ╠═161737a5-1f5b-4497-b00f-afc85da70569
# ╠═7252e0bc-c0cb-4477-a108-0896ae66c7ea
# ╠═6220c7db-be0c-47b6-b59e-748eabaadd87
# ╠═adaec625-d75e-429d-b6a2-18ec51c80c4a
# ╠═dd4b8dd3-f151-4564-af94-8102cbf6917f
# ╠═5604af2b-1435-4c64-92d8-546e758814bd
# ╠═99b96798-4a5a-4f75-a576-b03488487e02
# ╠═b756b05f-b4f5-4a64-b344-bf360e9a8bd3
# ╠═abe07005-c0b9-45d2-869c-680af950bf5f
# ╠═e4dc6b05-8f9d-4177-87a8-7c8ae622894b
# ╠═f7517c96-ae0d-43bb-9f12-f0d077172e70
# ╠═76bfd0d4-f723-44d8-b2a5-3f2645a9e8cb
# ╠═748afc4e-a8f7-4a12-99aa-296aff325520
# ╠═96b195c0-685e-48ec-a9f2-92e69e17b31c
# ╠═4f3ed2b1-4089-4caa-8180-d584b739d66f
# ╠═8c528821-ba3f-45e2-965d-e93aab90fc72
# ╠═b33cb49f-df54-44e9-af9e-bd2b6aa1e42a
# ╠═762b903f-0ded-4e3b-804f-59fc73c37dfa
# ╠═96349783-915d-4328-bc41-4e8d7e516ad7
# ╠═21966bd0-143a-4e95-9881-b645b9b5a261
# ╠═53653580-28c9-46ce-a9e8-63ec523e8022
# ╠═b4a38a4a-45ae-4666-803b-0769c0ea06de
# ╠═3ac42368-de6d-4012-98ab-cf275a9436a8
# ╠═be2b3c00-d909-48e9-b009-59f2fd996534
# ╠═dd238bab-06c5-453a-9a12-bd0711b82c42
# ╠═81926f04-dd75-4e09-b1b2-e39efbf48151
# ╠═111eca50-bfdd-42e7-b495-a4ee8cc6d48b
# ╠═d63552db-2eed-46e0-ad6a-801644e6ca2a
# ╠═ea18fef7-c1c6-46b9-b362-94e6762fc52f
# ╠═54ca7e3e-8b56-4b77-9488-f0c4ea98a1e4
# ╠═11478b83-6854-466c-8817-9fced3bc4b7b
# ╠═1ee2688f-0fb4-4227-b101-156b310a34ee
# ╠═49ca550c-5662-4ae3-841e-4411c670c6b7
# ╠═aab4cc7e-15c0-43a3-b9cb-35e1559478d1
# ╠═d82b6d02-a0ef-418b-898f-6dca245ff100
# ╠═d2b92466-6365-41b7-822e-2bf2df895dd3
# ╠═dcf44d38-95df-4e27-b4ef-3858cb0a40ac
# ╠═2c4015ad-8778-44f9-9d27-58a5ccbb5bd5
# ╠═b588e294-6449-45c7-afbc-6535f1dfa6eb
# ╠═e6028e8e-e58e-475e-994f-60bdecf164eb
# ╠═91597c20-dace-4a3b-854a-8682173d9476
# ╠═57fea3b4-eb6c-4064-87b5-314d83b1e117
# ╠═36eb9427-db56-4265-9621-0a91efe99a93

### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ c8ad2386-64d0-4318-a196-c3ea0bc62442
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 7feaa2b1-0dd8-4980-9b9d-e44b0991b6dc
using DataFrames: DataFrame

# ╔═╡ 778797cd-a008-4628-9ff1-e5aa947770d3
using BioSequences: LongRNA

# ╔═╡ 29b58ff3-7ea4-48bb-8f58-9859dff108b4
using Makie: @L_str

# ╔═╡ 94904ae5-cd2f-4e67-a2d7-8ac978b5ce15
using Statistics: mean, cor

# ╔═╡ a6a3644c-7307-4ddb-85d2-dc69e51adab6
using NaNStatistics: nanmean, nanstd

# ╔═╡ 352efc5c-45e3-11f0-0941-b591bfba1ed0
md"# Imports"

# ╔═╡ 9651cc83-074d-4c11-910e-93a392ce72c6
import Makie, CairoMakie, PlutoUI

# ╔═╡ d00f6685-e171-41d5-b887-bf12968d86eb
import FASTX, Infernal, Rfam, ViennaRNA

# ╔═╡ fd78654d-418e-420a-8eee-77498712057d
import SamApp2025

# ╔═╡ ba770877-fa6c-4bbe-ad6e-880ba37f133f
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ b1c0eadd-b87d-40fa-b2ac-5c56bec6533a
PlutoUI.TableOfContents()

# ╔═╡ 390fad57-e9fc-4833-83fa-2543bd675626
md"# Load data"

# ╔═╡ 1f271191-9cc1-4a52-9655-94eb5739d6bc
_hallmark_sites = SamApp2025.hallmark_sites_20230507

# ╔═╡ f40c1495-a22e-45d5-9494-caa4872374d6
ss_bps, ss_nps, ss_pks = SamApp2025.RF00162_sites_paired()

# ╔═╡ ccbb4d7f-4aa7-4d08-94ea-e45680c1c509
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20240730_with_pdb();

# ╔═╡ ebe2ba11-bac0-4796-ad91-b84d3cb29937
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 05e2322e-bfea-4138-b5f8-01a4478fa858
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 1c1dafbe-1583-49f3-b9d1-d7bf433180fa
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ a782b435-219e-4518-be19-5d653bc5a689
shape_rep0_aptamer_rbm_energies = [
    ismissing(seq) ? missing :
    RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq)))
    for seq = shape_data_045.aligned_sequences
];

# ╔═╡ a2a722dd-960c-4b71-b7c1-465b964691dc
shape_rep0_aptamer_rbm_energies[[299, 207]]

# ╔═╡ 9ce756f2-e1e2-49b3-9567-75ae6adaf369
shape_rep0_seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70");

# ╔═╡ cd6d9ad4-0a0d-498a-8505-eec983b6e4a1
shape_rep0_hits_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30");

# ╔═╡ 01ad33aa-b2b4-42f8-b1a3-c55d5b104073
shape_rep0_nats_seqs = shape_rep0_seed_seqs ∪ shape_rep0_hits_seqs;

# ╔═╡ 8999f767-0ee5-4bae-ad6f-726d07e09a8c
shape_rep0_rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm");

# ╔═╡ 5f37363c-bd93-4350-9e84-1218fd4ef219
shape_rep0_CM_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf");

# ╔═╡ a107fca3-82a6-496f-ba93-a71a0fe02b11
shape_rep0_rbmlo_seqs = shape_rep0_rbm_seqs ∩ findall((!ismissing).(shape_rep0_aptamer_rbm_energies) .&& (shape_rep0_aptamer_rbm_energies .< -300));

# ╔═╡ dba214be-a1d2-4f12-b26b-2f7ca5a046ce
ΔR_shape_rep0_sam = (
    nanmean(shape_data_rep0.shape_reactivities[:, :, conds_sam_rep0]; dim=3) .-
    shape_data_rep0.shape_reactivities[:, :, only(conds_mg_rep0)]
);

# ╔═╡ dfde59b3-7635-4cc9-b502-008abc31f87a
ΔR_shape_rep0_sam_avg_seed = nanmean(ΔR_shape_rep0_sam[:, shape_rep0_seed_seqs]; dim=2);

# ╔═╡ a9346418-f15c-444d-97e2-19b4f35ad2f1
ΔR_shape_rep0_sam_std_seed = nanstd(ΔR_shape_rep0_sam[:, shape_rep0_seed_seqs]; dim=2);

# ╔═╡ b7f4a268-0188-4051-83df-b171463b10e7
ΔR_shape_rep0_sam_avg_rbmlo = nanmean(ΔR_shape_rep0_sam[:, shape_rep0_rbmlo_seqs]; dim=2);

# ╔═╡ f7884b7e-e092-4bb6-9217-f5beae767d7b
ΔR_shape_rep0_sam_std_rbmlo = nanstd(ΔR_shape_rep0_sam[:, shape_rep0_rbmlo_seqs]; dim=2);

# ╔═╡ 85321c80-f7a5-46df-a648-61367243905a
ΔR_shape_rep0_sam_avg_CM = nanmean(ΔR_shape_rep0_sam[:, shape_rep0_CM_seqs]; dim=2)

# ╔═╡ 5e79a9e2-b5aa-4622-be90-e0521cee3022
ΔR_shape_rep0_sam_std_CM = nanstd(ΔR_shape_rep0_sam[:, shape_rep0_CM_seqs]; dim=2);

# ╔═╡ 8b67406b-501d-4a62-8cf9-2e27f19be1b0
shape_data_all_merged = SamApp2025.load_shapemapper_data_pierre_demux_20240801_with_pdb_repls_merged();

# ╔═╡ 50ddca2d-9734-4bfb-b170-3e804550b789
conds_SAM_all_merged = map(identity, indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_allrep", "SAMAP_1M7_1SAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 13a3c9b9-b1eb-496c-88c3-e46bfa22117f
conds_Mg_all_merged = map(identity, indexin(["SAMAP_1M7_noSAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 8b525135-2e4a-44f3-bade-16850ae91d65
shape_data_500 = SamApp2025.load_shapemapper_data_500v2_20240315();

# ╔═╡ 60f08a43-a24b-4c3d-a0f5-1cfa476b037c
aptamer_names_500 = "APSAM-S2-" .* lpad.(0:499, 3, '0')

# ╔═╡ 6a7cc1d5-4b97-427b-9e8b-745c3d4362bc
dms_df = SamApp2025.load_dms_data_sequences_table_20250303_with_aligned_sequences();

# ╔═╡ 635fd719-612e-4207-96d5-b4da5d415f9c
dms_data = SamApp2025.load_dms_data_20250303();

# ╔═╡ a3f04906-c7e6-49f4-8ee0-24530794671e
dms_aptamer_origin = [
	begin
		if n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_rbm"] ∪ aptamer_names_500[shape_data_500.aptamer_origin .== "rbm"]
			"rbm"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_seed70" .|| shape_data_045.aptamer_origin .== "RF00162_full30"]
			"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_inf"] ∪ aptamer_names_500[shape_data_500.aptamer_origin .== "infernal"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "Infrared"]
			"Infrared"
		end
	end for n = dms_data.aptamer_names
]

# ╔═╡ 93ae357f-ed40-4d07-bcf9-7f86c6274e9e
dms_aptamer_rbm_energies = [ismissing(seq) ? missing : RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq))) for seq = dms_data.aligned_sequence];

# ╔═╡ c099965d-e0d3-48c1-87c1-b43cf488455f
dms_nat_seqs = findall(dms_aptamer_origin .== "natural")

# ╔═╡ 71ac67b0-1167-4256-8fb3-fe4490314e81
dms_rbmlo_seqs = findall(dms_aptamer_origin .== "rbm") ∩ findall((!ismissing).(dms_aptamer_rbm_energies) .&& (dms_aptamer_rbm_energies .< -300))

# ╔═╡ e8f4d853-bfd4-4976-be85-3a74d0aa5185
dms_rbmhi_seqs = findall(dms_aptamer_origin .== "rbm") ∩ findall((!ismissing).(dms_aptamer_rbm_energies) .&& (dms_aptamer_rbm_energies .> -300))

# ╔═╡ fddb0d6b-89aa-4769-a3cb-4e002064e858
conds_sam_dms = [1]

# ╔═╡ f986a816-a856-4953-b9b7-e496726600f2
conds_mg_dms = [2]

# ╔═╡ c388c632-c9cb-4069-8c59-f680ce1aab15
ΔR_dms_sam = dms_data.shape_reactivities[:, :, only(conds_sam_dms)] - dms_data.shape_reactivities[:, :, only(conds_mg_dms)];

# ╔═╡ a8bade84-8af2-4d01-9a7f-1600ccceb37d
ΔR_dms_sam_avg_natural = [nanmean([ΔR_dms_sam[i,n] for n = dms_nat_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 7eda1c34-3f5e-4a2a-8fbc-df84962b2b57
ΔR_dms_sam_std_natural = [nanstd([ΔR_dms_sam[i,n] for n = dms_nat_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 077cfd71-12d2-480e-96f2-94904fc02719
ΔR_dms_sam_avg_rbmlo = [nanmean([ΔR_dms_sam[i,n] for n = dms_rbmlo_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 1183eeb9-6cdb-467d-b663-ac4748495ee1
ΔR_dms_sam_std_rbmlo = [nanstd([ΔR_dms_sam[i,n] for n = dms_rbmlo_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 298593d9-2576-47dd-9a1d-69f10379a1c3
ΔR_dms_sam_avg_rbmhi = [nanmean([ΔR_dms_sam[i,n] for n = dms_rbmhi_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 614620b9-6dc2-477f-992c-4a42682bb68b
ΔR_dms_sam_std_rbmhi = [nanstd([ΔR_dms_sam[i,n] for n = dms_rbmhi_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ b77793f8-2495-40cf-a09a-5c39defdac1c
dms_aptamers_AC_content = [mean(seq[i] ∈ ('A', 'C') for seq = skipmissing(dms_data.aligned_sequence)) for i = 1:108]

# ╔═╡ 40c952e1-20cd-4e1a-81bf-374b5197b291
md"# Preliminaries"

# ╔═╡ 59ab1fa3-95c4-475f-95c0-7be5d93d1a1f
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

# ╔═╡ c4b6b4d2-d19f-4c31-93f3-b55c7339c896
md"# Figure 5"

# ╔═╡ cb7081bf-6a40-44d3-a24a-66f5b5474926
let fig = Makie.Figure()

	width = 900
	profiles_height = 85
	avg_delta_height = 120
	profiles_ylabel = "react."
	profiles_Δ_ylabel = "Δreact."

	function dummy_ax(row)
		_dummy_ax = Makie.Axis(fig[row,1]; height=10, width=20, xgridvisible=false, ygridvisible=false)
		Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
		Makie.hidexdecorations!(_dummy_ax)
		Makie.hideydecorations!(_dummy_ax)
		return _dummy_ax
	end

	xticks = 5:5:108

	natural_examples = [
		only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), # 4KQY
		only(findall(shape_data_045.aptamer_names .== "APSAMN7"))
	]

	@show shape_data_045.aptamer_ids[natural_examples[1]] shape_data_045.aptamer_ids[natural_examples[2]]

	axes_shape_reactivity_profile_natural_examples = [
		Makie.Axis(fig[row,1]; valign=:bottom, width, height=profiles_height, xticks, ylabel=profiles_ylabel, xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true) for row=[1,4]
	]
	axes_shape_reactivity_diff_natural_examples = [
		Makie.Axis(fig[row,1]; valign=:bottom, width, height=profiles_height, xticks, xlabel="site", ylabel=profiles_Δ_ylabel, xgridvisible=false, ygridvisible=false, yticks=[-2,-1,0], xtrimspine=true, ytrimspine=true) for row = [2, 5]
	]
	# add some space in between
	dummy_ax(3)

	for (n_ex, ax_R, ax_Δ) = zip(natural_examples, axes_shape_reactivity_profile_natural_examples, axes_shape_reactivity_diff_natural_examples)
		shape_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
		shape_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
		for (x0, xf, color, alpha) = struct_bands, ax = (ax_R, ax_Δ)
		    Makie.vspan!(ax, x0, xf; color=(color, alpha))
		end
		Makie.stairs!(ax_R, 1:108, shape_R_mg; step=:center, color=:gray, label="no SAM")
		Makie.stairs!(ax_R, 1:108, shape_R_sam; step=:center, color=:purple, label="with SAM")
		Makie.barplot!(ax_Δ, 1:108, shape_R_sam - shape_R_mg; color=ifelse.(shape_R_sam - shape_R_mg .< 0, :green, :gray))
		#Makie.scatter!(ax_Δ, SamApp2025.hallmark_sites_20230507, -1.9one.(SamApp2025.hallmark_sites_20230507), markersize=7, color=:black, marker=:utriangle)
		Makie.scatter!(ax_Δ, SamApp2025.hallmark_sites_20230507, fill(minimum([x for x = shape_R_sam - shape_R_mg if !isnan(x)]), length(SamApp2025.hallmark_sites_20230507)); markersize=7, color=:black, marker=:utriangle)
		for ax = (ax_Δ, ax_R)
			Makie.hidespines!(ax, :r, :t, :b)
			Makie.hidexdecorations!(ax)
			Makie.xlims!(ax, 0.5, 108.5)
		end
		# Makie.ylims!(ax_Δ, -2, 0.5)
		# Makie.ylims!(ax_R, -0.5, 6)
	end

	_xs = 1:108

	_std_thresh = 0
	ΔR_dms_sam_avg_natural_no0std = [ΔR_dms_sam_std_natural[i] > _std_thresh ? ΔR_dms_sam_avg_natural[i] : NaN for i = eachindex(ΔR_dms_sam_avg_natural)]
	ΔR_dms_sam_avg_rbmlo_no0std = [ΔR_dms_sam_std_rbmlo[i] > _std_thresh ? ΔR_dms_sam_avg_rbmlo[i] : NaN for i = eachindex(ΔR_dms_sam_std_rbmlo)]
	ΔR_dms_sam_avg_rbmhi_no0std = [ΔR_dms_sam_std_rbmhi[i] > _std_thresh ? ΔR_dms_sam_avg_rbmhi[i] : NaN for i = eachindex(ΔR_dms_sam_avg_rbmhi)]

	ax_shape_avg_reactivity_profile_natural_only = Makie.Axis(fig[6,:]; width, height=avg_delta_height, xticks=5:10:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, ylabel="⟨Δr⟩ (SHAPE)", xtrimspine=true, ytrimspine=true)
	Makie.band!(ax_shape_avg_reactivity_profile_natural_only, _xs, (ΔR_shape_rep0_sam_avg_seed - ΔR_shape_rep0_sam_std_seed/2)[_xs], (ΔR_shape_rep0_sam_avg_seed + ΔR_shape_rep0_sam_std_seed/2)[_xs], color=(:gray, 0.25))
	Makie.scatterlines!(ax_shape_avg_reactivity_profile_natural_only, _xs, ΔR_shape_rep0_sam_avg_seed[_xs], markersize=5, linewidth=0.5, color=:black, label="Natural")
	Makie.xlims!(ax_shape_avg_reactivity_profile_natural_only, 0.5, 108.5)
	Makie.hidespines!(ax_shape_avg_reactivity_profile_natural_only, :t, :r, :b)
	Makie.hidexdecorations!(ax_shape_avg_reactivity_profile_natural_only)

	ax_dms_avg_reactivity_profile_natural_only = Makie.Axis(fig[7,:]; width, height=avg_delta_height, xticks=5:10:108, yticks=-2:0.5:1, xgridvisible=false, ygridvisible=false, ylabel="⟨Δr⟩ (DMS)", xtrimspine=true, ytrimspine=true)
	Makie.band!(ax_dms_avg_reactivity_profile_natural_only, _xs, (ΔR_dms_sam_avg_natural_no0std - ΔR_dms_sam_std_natural/2)[_xs], (ΔR_dms_sam_avg_natural_no0std + ΔR_dms_sam_std_natural/2)[_xs]; color=(:gray, 0.25))
	Makie.scatterlines!(ax_dms_avg_reactivity_profile_natural_only, _xs, ΔR_dms_sam_avg_natural_no0std; linewidth=0.5, markersize=5, color=:black)
	#Makie.axislegend(ax_dms_avg_reactivity_profile_natural_only, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax_dms_avg_reactivity_profile_natural_only, :t, :r, :b)
	Makie.hidexdecorations!(ax_dms_avg_reactivity_profile_natural_only)
	Makie.xlims!(ax_dms_avg_reactivity_profile_natural_only, 0.5, 108.5)
	Makie.ylims!(ax_dms_avg_reactivity_profile_natural_only, -0.7, 0.7)

	ax_AC = Makie.Axis(fig[8,:]; width, height=40, xticks, yticks=[0,1], xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="A,C freq.", xtrimspine=true, ytrimspine=true)
	Makie.barplot!(ax_AC, dms_aptamers_AC_content; width=0.5, color=:lightblue)
	Makie.xlims!(ax_AC, 0.5, 108.5)
	Makie.hidespines!(ax_AC, :t, :r)
	#Makie.hidexdecorations!(ax_AC)
	dummy_ax(9)

	Makie.linkxaxes!(
		axes_shape_reactivity_profile_natural_examples..., axes_shape_reactivity_diff_natural_examples...,
		ax_dms_avg_reactivity_profile_natural_only, ax_AC
	)
	# Makie.Label(fig[1,1][1,1,Makie.TopLeft()], "A)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,2][1,1,Makie.TopLeft()], "B)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,3][1,1,Makie.TopLeft()], "C)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[2,:][1,1,Makie.TopLeft()], "D)", font=:bold, padding=(0,0,0,0))
	# Makie.Label(fig[3,:][1,1,Makie.TopLeft()], "E)", font=:bold, padding=(0,0,0,0))

	Makie.resize_to_layout!(fig)
	save_path = mktempdir() * "/fig.pdf"
	println("Saving to $save_path")
	Makie.save(save_path, fig)
	fig
end

# ╔═╡ fc78947b-3c85-4414-a28e-ac3f9401a541
only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")) # 4KQY

# ╔═╡ 330cb974-90cf-4b0e-ad5b-2b7480199b92
shape_data_045.aligned_sequences[308]

# ╔═╡ db8dd686-4375-44b2-8dbb-7b6664932a18
99 ∈ SamApp2025.hallmark_sites_20230507

# ╔═╡ e3300392-64b4-439c-9373-cbbf82e7083b
md"# Figure 6"

# ╔═╡ 2a83ca52-1a19-47e1-b7fb-522e333193bd
let fig = Makie.Figure()

	width = 900
	profiles_height = 85
	avg_delta_height = 120
	profiles_ylabel = "react."
	profiles_Δ_ylabel = "Δreact."

	function dummy_ax(row)
		_dummy_ax = Makie.Axis(fig[row,1]; height=10, width=20, xgridvisible=false, ygridvisible=false)
		Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
		Makie.hidexdecorations!(_dummy_ax)
		Makie.hideydecorations!(_dummy_ax)
		return _dummy_ax
	end

	xticks = 5:5:108

	RBM_examples = 	[299, 207]

	_xs = 1:108

	_std_thresh = 0
	ΔR_dms_sam_avg_natural_no0std = [ΔR_dms_sam_std_natural[i] > _std_thresh ? ΔR_dms_sam_avg_natural[i] : NaN for i = eachindex(ΔR_dms_sam_avg_natural)]
	ΔR_dms_sam_avg_rbmlo_no0std = [ΔR_dms_sam_std_rbmlo[i] > _std_thresh ? ΔR_dms_sam_avg_rbmlo[i] : NaN for i = eachindex(ΔR_dms_sam_std_rbmlo)]
	ΔR_dms_sam_avg_rbmhi_no0std = [ΔR_dms_sam_std_rbmhi[i] > _std_thresh ? ΔR_dms_sam_avg_rbmhi[i] : NaN for i = eachindex(ΔR_dms_sam_avg_rbmhi)]

	axes_shape_reactivity_profile_RBM_examples = [
		Makie.Axis(fig[row,1]; width, height=profiles_height, xticks, ylabel="reactivity", xgridvisible=false, ygridvisible=false, yticks=0:2:5, ytrimspine=true) for row = [1,4]
	]
	axes_shape_reactivity_diff_RBM_examples = [
		Makie.Axis(fig[row,1]; width, height=profiles_height, xticks, xlabel="site", ylabel="Δreactivity", xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true, yticks=-2:2) for row=[2,5]
	]
	dummy_ax(3) # add some space in between
	for (n_ex, ax_R, ax_Δ) = zip(RBM_examples, axes_shape_reactivity_profile_RBM_examples, axes_shape_reactivity_diff_RBM_examples)
		for (x0, xf, color, alpha) = struct_bands, ax = (ax_R, ax_Δ)
		    Makie.vspan!(ax, x0, xf; color=(color, alpha))
		end
		_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
		_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
		Makie.stairs!(ax_R, 1:108, _R_mg, color=:gray, step=:center)
		Makie.stairs!(ax_R, 1:108, _R_sam, color=:purple, step=:center)
		Makie.barplot!(ax_Δ, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
		Makie.scatter!(ax_Δ, SamApp2025.hallmark_sites_20230507, fill(minimum([x for x = _R_sam - _R_mg if !isnan(x)]), length(SamApp2025.hallmark_sites_20230507)), markersize=7, color=:black, marker=:utriangle)
		# Makie.ylims!(ax_Δ, -1.5, 1.5)
		# Makie.ylims!(ax_R, -1, 4)
		for ax = (ax_R, ax_Δ)
			Makie.xlims!(ax, 0.5, 108.5)
			Makie.hidespines!(ax, :t, :r, :b)
			Makie.hidexdecorations!(ax)
		end
	end

	_xs = 3:107

	ax_shape_avg_RBM = Makie.Axis(fig[6,:]; width, height=avg_delta_height, xticks=5:5:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, ylabel="⟨Δreactivity⟩", xtrimspine=true, ytrimspine=true)
	Makie.band!(ax_shape_avg_RBM, _xs, (ΔR_shape_rep0_sam_avg_seed - ΔR_shape_rep0_sam_std_seed/2)[_xs], (ΔR_shape_rep0_sam_avg_seed + ΔR_shape_rep0_sam_std_seed/2)[_xs], color=(:gray, 0.25))
	Makie.band!(ax_shape_avg_RBM, _xs, (ΔR_shape_rep0_sam_avg_rbmlo - ΔR_shape_rep0_sam_std_rbmlo/2)[_xs], (ΔR_shape_rep0_sam_avg_rbmlo + ΔR_shape_rep0_sam_std_rbmlo/2)[_xs], color=(:blue, 0.25))
	Makie.scatterlines!(ax_shape_avg_RBM, _xs, ΔR_shape_rep0_sam_avg_seed[_xs]; markersize=5, linewidth=1, color=:gray)
	Makie.scatterlines!(ax_shape_avg_RBM, _xs, ΔR_shape_rep0_sam_avg_rbmlo[_xs]; markersize=5, linewidth=1, color=:blue)
	#Makie.axislegend(ax_shape_avg_RBM, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.xlims!(0.5, 108.5)
	Makie.hidespines!(ax_shape_avg_RBM, :t, :r, :b)
	Makie.hidexdecorations!(ax_shape_avg_RBM)

	ax_shape_avg_CM = Makie.Axis(fig[7,:]; width, height=avg_delta_height, xticks=5:5:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="⟨Δreactivity⟩", xtrimspine=true, ytrimspine=true)
	Makie.band!(ax_shape_avg_CM, _xs, (ΔR_shape_rep0_sam_avg_seed - ΔR_shape_rep0_sam_std_seed/2)[_xs], (ΔR_shape_rep0_sam_avg_seed + ΔR_shape_rep0_sam_std_seed/2)[_xs], color=(:gray, 0.25))
	Makie.band!(ax_shape_avg_CM, _xs, (ΔR_shape_rep0_sam_avg_CM - ΔR_shape_rep0_sam_std_CM/2)[_xs], (ΔR_shape_rep0_sam_avg_CM + ΔR_shape_rep0_sam_std_CM/2)[_xs], color=(:red, 0.25))
	Makie.scatterlines!(ax_shape_avg_CM, _xs, ΔR_shape_rep0_sam_avg_seed[_xs]; markersize=5, linewidth=1, color=:gray)
	Makie.scatterlines!(ax_shape_avg_CM, _xs, ΔR_shape_rep0_sam_avg_CM[_xs]; markersize=5, linewidth=1, color=:red)
	#Makie.scatter!(ax_shape_avg_CM, SamApp2025.hallmark_sites_20230507, -2.0one.(SamApp2025.hallmark_sites_20230507), markersize=7, color=:black, marker=:utriangle)
	#Makie.axislegend(ax_shape_avg_CM, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax_shape_avg_CM, :t, :r)
	Makie.xlims!(0.5, 108.5)

	Makie.linkxaxes!(
		axes_shape_reactivity_profile_RBM_examples..., axes_shape_reactivity_diff_RBM_examples...,
		ax_shape_avg_RBM, ax_shape_avg_CM
	)
	# Makie.Label(fig[1,1][1,1,Makie.TopLeft()], "A)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,2][1,1,Makie.TopLeft()], "B)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,3][1,1,Makie.TopLeft()], "C)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[2,:][1,1,Makie.TopLeft()], "D)", font=:bold, padding=(0,0,0,0))
	# Makie.Label(fig[3,:][1,1,Makie.TopLeft()], "E)", font=:bold, padding=(0,0,0,0))

	Makie.resize_to_layout!(fig)
	save_path = mktempdir() * "/fig.pdf"
	println("Saving to $save_path")
	Makie.save(save_path, fig)
	fig
end

# ╔═╡ c386b292-ed69-4231-a1fb-a9f8b4f86804
md"# Figure 7"

# ╔═╡ b27a2a43-3eae-443b-9e4e-a7fa18db47bc
let fig = Makie.Figure()
	hist_size = 200

	ax = Makie.Axis(
		fig[1,1], width=hist_size, height=hist_size, xlabel="SHAPE reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true
	)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_bps, shape_rep0_nats_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_nps, shape_rep0_nats_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_bps, shape_rep0_nats_seqs, conds_sam_rep0])), label="base paired", normalization=:pdf, bins=-2:0.05:6, linewidth=2, color=:teal)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_nps, shape_rep0_nats_seqs, conds_sam_rep0])), label="not paired", normalization=:pdf, bins=-2:0.05:6, linewidth=2, color=:orange)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	#Makie.axislegend(ax, framevisible=false, patchlabelgap=3, position=(-0.02, 1))
	#Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r)

	_dummy_ax = Makie.Axis(fig[1,2], width=10, xgridvisible=false, ygridvisible=false)
	Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
	Makie.hidexdecorations!(_dummy_ax)
	Makie.hideydecorations!(_dummy_ax)

	ax = Makie.Axis(fig[1,3], width=hist_size, height=hist_size, xlabel="SHAPE reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_bps, shape_rep0_nats_seqs, conds_sam_rep0])), label="b.p.", normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_nps, shape_rep0_nats_seqs, conds_sam_rep0])), label="n.p.", normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_pks, shape_rep0_nats_seqs, conds_mg_rep0])), label="p.k.", normalization=:pdf, bins=-2:0.1:6, linewidth=2, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)

	ax = Makie.Axis(fig[1,4], width=hist_size, height=hist_size, xlabel="SHAPE reactivity", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_bps, shape_rep0_nats_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_nps, shape_rep0_nats_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[ss_pks, shape_rep0_nats_seqs, conds_sam_rep0])), label="pseudoknot", normalization=:pdf, bins=-2:0.1:6, linewidth=2, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	#Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)


	# DMS

	bp_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = dms_nat_seqs for i = ss_bps for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
	np_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = dms_nat_seqs for i = ss_nps for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

	pk_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = dms_nat_seqs for i = ss_pks for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
	pk_nat_mg_reactivities = [dms_data.shape_reactivities[i,n,c] for n = dms_nat_seqs for i = ss_pks for c = conds_mg_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

	ax = Makie.Axis(fig[2,1], width=hist_size, height=hist_size, xlabel="DMS reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:2:10, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), normalization=:pdf, bins=-1:0.05:2, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), normalization=:pdf, bins=-1:0.05:2, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), label="base paired", normalization=:pdf, bins=-1:0.05:2, linewidth=2, color=:teal)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), label="not paired", normalization=:pdf, bins=-1:0.05:2, linewidth=2, color=:orange)
	Makie.xlims!(-0.5, 2)
	Makie.ylims!(-0.07, 10)
	#Makie.axislegend(ax, framevisible=false, patchlabelgap=3, position=(-0.02, 1))
	#Makie.axislegend(ax, position=(1, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r)

	_dummy_ax = Makie.Axis(fig[2,2], width=10, xgridvisible=false, ygridvisible=false)
	Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
	Makie.hidexdecorations!(_dummy_ax)
	Makie.hideydecorations!(_dummy_ax)

	ax = Makie.Axis(fig[2,3], width=hist_size, height=hist_size, xlabel="DMS reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), label="b.p.", normalization=:pdf, bins=-1:0.05:2, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), label="n.p.", normalization=:pdf, bins=-1:0.05:2, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, pk_nat_mg_reactivities), label="pseudoknot", normalization=:pdf, bins=-1:0.05:2, linewidth=2, color=:black)
	Makie.xlims!(-0.5, 2)
	Makie.ylims!(-0.07, 10)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)

	ax = Makie.Axis(fig[2,4], width=hist_size, height=hist_size, xlabel="DMS reactivity", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, pk_nat_sam_reactivities), label="pseudoknot", normalization=:pdf, bins=-1:0.05:2, linewidth=2, color=:black)
	Makie.xlims!(-0.5, 2)
	Makie.ylims!(-0.07, 10)
	#Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)

	Makie.resize_to_layout!(fig)
	save_path = mktempdir() * "/fig.pdf"
	println("Saving to $save_path")
	Makie.save(save_path, fig)
	fig
end

# ╔═╡ Cell order:
# ╠═352efc5c-45e3-11f0-0941-b591bfba1ed0
# ╠═c8ad2386-64d0-4318-a196-c3ea0bc62442
# ╠═9651cc83-074d-4c11-910e-93a392ce72c6
# ╠═d00f6685-e171-41d5-b887-bf12968d86eb
# ╠═fd78654d-418e-420a-8eee-77498712057d
# ╠═ba770877-fa6c-4bbe-ad6e-880ba37f133f
# ╠═7feaa2b1-0dd8-4980-9b9d-e44b0991b6dc
# ╠═778797cd-a008-4628-9ff1-e5aa947770d3
# ╠═29b58ff3-7ea4-48bb-8f58-9859dff108b4
# ╠═94904ae5-cd2f-4e67-a2d7-8ac978b5ce15
# ╠═a6a3644c-7307-4ddb-85d2-dc69e51adab6
# ╠═b1c0eadd-b87d-40fa-b2ac-5c56bec6533a
# ╠═390fad57-e9fc-4833-83fa-2543bd675626
# ╠═1f271191-9cc1-4a52-9655-94eb5739d6bc
# ╠═f40c1495-a22e-45d5-9494-caa4872374d6
# ╠═ccbb4d7f-4aa7-4d08-94ea-e45680c1c509
# ╠═ebe2ba11-bac0-4796-ad91-b84d3cb29937
# ╠═05e2322e-bfea-4138-b5f8-01a4478fa858
# ╠═1c1dafbe-1583-49f3-b9d1-d7bf433180fa
# ╠═a782b435-219e-4518-be19-5d653bc5a689
# ╠═a2a722dd-960c-4b71-b7c1-465b964691dc
# ╠═9ce756f2-e1e2-49b3-9567-75ae6adaf369
# ╠═cd6d9ad4-0a0d-498a-8505-eec983b6e4a1
# ╠═01ad33aa-b2b4-42f8-b1a3-c55d5b104073
# ╠═8999f767-0ee5-4bae-ad6f-726d07e09a8c
# ╠═5f37363c-bd93-4350-9e84-1218fd4ef219
# ╠═a107fca3-82a6-496f-ba93-a71a0fe02b11
# ╠═dba214be-a1d2-4f12-b26b-2f7ca5a046ce
# ╠═dfde59b3-7635-4cc9-b502-008abc31f87a
# ╠═a9346418-f15c-444d-97e2-19b4f35ad2f1
# ╠═b7f4a268-0188-4051-83df-b171463b10e7
# ╠═f7884b7e-e092-4bb6-9217-f5beae767d7b
# ╠═85321c80-f7a5-46df-a648-61367243905a
# ╠═5e79a9e2-b5aa-4622-be90-e0521cee3022
# ╠═8b67406b-501d-4a62-8cf9-2e27f19be1b0
# ╠═50ddca2d-9734-4bfb-b170-3e804550b789
# ╠═13a3c9b9-b1eb-496c-88c3-e46bfa22117f
# ╠═8b525135-2e4a-44f3-bade-16850ae91d65
# ╠═60f08a43-a24b-4c3d-a0f5-1cfa476b037c
# ╠═6a7cc1d5-4b97-427b-9e8b-745c3d4362bc
# ╠═635fd719-612e-4207-96d5-b4da5d415f9c
# ╠═a3f04906-c7e6-49f4-8ee0-24530794671e
# ╠═93ae357f-ed40-4d07-bcf9-7f86c6274e9e
# ╠═c099965d-e0d3-48c1-87c1-b43cf488455f
# ╠═71ac67b0-1167-4256-8fb3-fe4490314e81
# ╠═e8f4d853-bfd4-4976-be85-3a74d0aa5185
# ╠═fddb0d6b-89aa-4769-a3cb-4e002064e858
# ╠═f986a816-a856-4953-b9b7-e496726600f2
# ╠═c388c632-c9cb-4069-8c59-f680ce1aab15
# ╠═a8bade84-8af2-4d01-9a7f-1600ccceb37d
# ╠═7eda1c34-3f5e-4a2a-8fbc-df84962b2b57
# ╠═077cfd71-12d2-480e-96f2-94904fc02719
# ╠═1183eeb9-6cdb-467d-b663-ac4748495ee1
# ╠═298593d9-2576-47dd-9a1d-69f10379a1c3
# ╠═614620b9-6dc2-477f-992c-4a42682bb68b
# ╠═b77793f8-2495-40cf-a09a-5c39defdac1c
# ╠═40c952e1-20cd-4e1a-81bf-374b5197b291
# ╠═59ab1fa3-95c4-475f-95c0-7be5d93d1a1f
# ╠═c4b6b4d2-d19f-4c31-93f3-b55c7339c896
# ╠═cb7081bf-6a40-44d3-a24a-66f5b5474926
# ╠═fc78947b-3c85-4414-a28e-ac3f9401a541
# ╠═330cb974-90cf-4b0e-ad5b-2b7480199b92
# ╠═db8dd686-4375-44b2-8dbb-7b6664932a18
# ╠═e3300392-64b4-439c-9373-cbbf82e7083b
# ╠═2a83ca52-1a19-47e1-b7fb-522e333193bd
# ╠═c386b292-ed69-4231-a1fb-a9f8b4f86804
# ╠═b27a2a43-3eae-443b-9e4e-a7fa18db47bc

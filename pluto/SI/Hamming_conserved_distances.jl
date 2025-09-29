### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 1dec11ae-386c-11ef-1b4c-51f221040ed7
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ eb41dfa5-12cc-4262-b0ce-15ed56c3a2cb
using NaNStatistics: nansum

# ╔═╡ e9a45fb3-7b6f-4957-8f8b-240c964e3fb6
using BioSequences: LongRNA

# ╔═╡ 072312ba-58f2-4bea-898c-5164c976ea16
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 0e9405c3-452c-4278-8b31-8aa89af09d06
using Statistics: mean

# ╔═╡ 8d4bc108-2cb6-416c-9a9d-e5cdf2120206
using LogExpFunctions: xlogx

# ╔═╡ b260395b-0f33-4289-83d8-8fdfa851c78c
md"""
# Imports
"""

# ╔═╡ d43e0635-c1f2-48ff-b92f-ae7000f305a3
import Makie, CairoMakie

# ╔═╡ 1786306b-e9a5-4013-a1f6-a1133c3f99e3
import Infernal, Rfam

# ╔═╡ 96134a6b-a910-4831-9c75-e7d052b04481
import SamApp2025

# ╔═╡ 8798da86-73ba-4474-81dd-10296384e48e
@show Rfam.get_rfam_directory() Rfam.get_rfam_version();

# ╔═╡ 48a3772d-1f32-490c-9978-3099beb188dc
md"""# Functions"""

# ╔═╡ 9853b581-ca77-4b4d-ac8b-6ca3a960d244
xlog2x(x) = xlogx(x) / log(oftype(x,2))

# ╔═╡ 3cc2bb41-3213-42d4-995e-736eeaf01534
md"""# Load data (Repl. 0)"""

# ╔═╡ 26683581-756c-47c3-98bf-91a5035875d0
# load SHAPE data
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 4d577148-6793-493c-be1e-50d447224f47
# split rep0 from rep4+5
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 2f5a46ec-7ab2-444e-b050-b2c2c45a7e4d
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 8d4140bf-d8cf-45fa-8ced-4bd7900a1721
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 6e530ef6-e22c-4d9d-a7d3-523682526e09
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 2326a0b9-2373-4629-b67e-f190f6d37102
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 9483c0e6-2fa1-489e-bfdc-625b766bc2a9
(; bps, nps, pks) = SamApp2025.RF00162_sites_paired()

# ╔═╡ 4ab7bcfc-5d88-4820-a623-e59f27550cde
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 58f69627-684e-4951-a447-57240bb9e38a
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ eaec262e-d68c-49f5-b14a-d4c1e6a10ee1
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 3af22d1d-dc5b-45e9-953f-054320b33deb
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ b5bd2ee2-289c-4ae8-8165-7ea81fb2bdfb
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ fb2046e6-a7fb-4ebf-9b0b-d34865c5f9a9
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ fbcd5369-f492-4857-bb30-7fa695bd06c5
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ c43cf4ff-1884-427b-9ac0-d4f9efa75a8a
shape_stats_rep0 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 33775511-874f-415a-8486-1bb317973ad0
_thresh = log(5)

# ╔═╡ 1b3278e0-231c-47de-ad66-068551a084c7
_sites = SamApp2025.hallmark_sites_20230507;

# ╔═╡ 609d1eea-270b-4ffd-8706-9b14849cc76a
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3));

# ╔═╡ efa6e6c0-f018-4c4c-9bca-8267947265f1
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3));

# ╔═╡ a4f602b8-5a0b-4147-9643-5e81aba14f48
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 4f39ca7b-75ab-4168-9314-64a06d1779e7
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ 376848db-5682-47d2-8279-0eb8b30e1bce
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ 1f2365e4-f6ef-449a-b244-073b4099fb3b
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

# ╔═╡ cafde511-5ef8-4684-a945-90139bb20826
aptamer_rbm_energies_rep0 = [
    ismissing(seq) ? missing :
    free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ 94439257-4a29-4a80-bb1c-e9afc7da3753
md"""# Load data (500)"""

# ╔═╡ 612e5568-84c5-4b8b-a94d-3277c05eac0d
shape_data_500 = SamApp2025.load_shapemapper_data_500v2_20240315();

# ╔═╡ 5c66a219-fabf-4187-aad0-fd7ddcdd675b
conds_sam_500 = [1,2];

# ╔═╡ a6adc734-4075-45b1-b11b-eee8b5d91677
conds_mg_500 = [4];

# ╔═╡ de88df39-55c1-4518-8a9e-82e66b471a3d
conds_30C_500 = [6];

# ╔═╡ d28b51e3-c24d-4636-958c-cdf66fccd09b
aptamer_rbm_energies_500 = free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(shape_data_500.aligned_sequences));

# ╔═╡ 120eac94-3480-4ba7-ba43-15a9fb84856c
shape_stats_500 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = shape_data_500.shape_reactivities[bps, :, conds_sam_500],
    unpaired_reactivities = shape_data_500.shape_reactivities[nps, :, conds_sam_500],
    all_reactivities = shape_data_500.shape_reactivities[:, :, conds_sam_500],
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ da8d155e-bd5b-4eff-a830-378474d0d9c4
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_mg_500]; dim=(1,3));

# ╔═╡ 88a27614-0463-44ee-8508-76b65e0936b9
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ 60ae1097-6576-40d6-a813-a2ff4f7975ec
_responds_sam_yes_500 = (x_mg_500 .< -_thresh) .& (x_sam_500 .> +_thresh);

# ╔═╡ 7ce05caa-e694-4a61-8085-7ffc3a732667
_responds_sam_nop_500 = (x_mg_500 .> +_thresh) .| (x_sam_500 .< -_thresh);

# ╔═╡ 1f090e07-cd4e-4ded-8da9-0d3be10451c2
_inconclusive_500 = ((!).(_responds_sam_yes_500)) .& ((!).(_responds_sam_nop_500));

# ╔═╡ 8a60ae44-7aad-4265-93c1-1caba00afd13
rbm_seqs_500 = findall(shape_data_500.aptamer_origin .== "rbm");

# ╔═╡ 7dcb4fff-c8ab-4f3a-b36c-19b9db528142
md"""
# Compute distances to natural sequences
"""

# ╔═╡ 8ec5db2f-da12-48f0-8e1e-7a8c97901742
RF00162_hits = SamApp2025.rfam_RF00162_hits();

# ╔═╡ 2d45af78-aa2e-4d8f-8fea-1e8e865abcb6
aptamer_natural_distances_500 = SamApp2025.hamming(SamApp2025.onehot(shape_data_500.aligned_sequences), SamApp2025.onehot(RF00162_hits));

# ╔═╡ bdec681b-179f-4e7e-9b5d-b2a4ca01d5af
aptamer_natural_distances_rep0 = [ismissing(seq) ? missing : minimum(SamApp2025.hamming(SamApp2025.onehot(LongRNA{4}(seq)), SamApp2025.onehot(RF00162_hits))) for seq = shape_data_045.aligned_sequences];

# ╔═╡ 87f6e1da-8490-45b3-8eff-6a78bc6a7bc9
p_site_hits = reshape(mean(SamApp2025.onehot(RF00162_hits); dims=3), 5, 108)

# ╔═╡ 5abad490-b3d9-4fb2-a87c-89b160c32b0b
poccupancy = p_site_hits ./ sum(p_site_hits; dims=1)

# ╔═╡ 01af5f88-a0ed-4ee5-8c8a-b744463f4313
conservation = vec((log2(5) .+ sum(xlog2x.(poccupancy); dims=1)));

# ╔═╡ b17a6262-e87e-4344-aa41-35cf78770d86
conserved_sites = conservation .> 0.5

# ╔═╡ fbac32fe-ecc7-443b-8f2a-6e0935531cba
sum(conserved_sites)

# ╔═╡ 09d81eb5-4759-486c-8e0a-71f7f56fb7cb
aptamer_natural_distances_500_conserved = SamApp2025.hamming(SamApp2025.onehot(shape_data_500.aligned_sequences)[:, conserved_sites, :], SamApp2025.onehot(RF00162_hits)[:, conserved_sites, :]);

# ╔═╡ f1700a78-538c-41d1-bb60-d644e0b27cec
aptamer_natural_distances_rep0_conserved = [ismissing(seq) ? missing : minimum(SamApp2025.hamming(SamApp2025.onehot(LongRNA{4}(seq))[:, conserved_sites, :], SamApp2025.onehot(RF00162_hits)[:, conserved_sites, :])) for seq = shape_data_045.aligned_sequences];

# ╔═╡ b9689ab6-a819-4f22-b8bc-1add6e7366c7
md"""# Plot"""

# ╔═╡ 7aa22d92-c726-401f-ad37-e4efd683df6e
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1], width=300, height=300, xlabel="Divergence from closest natural", ylabel="RBM score", title="Total distance", xticks=0:0.1:0.6, yticks=200:25:350, xgridvisible=false, ygridvisible=false)

	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500[rbm_seqs_500 ∩ findall(_inconclusive_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_inconclusive_500)]; markersize=15, color=(:gray, 0.5), label="RBM (inconcl.)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500[rbm_seqs_500 ∩ findall(_responds_sam_yes_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_responds_sam_yes_500)]; markersize=10, color=:blue, marker='●', label="RBM (✓)")
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500[rbm_seqs_500 ∩ findall(_responds_sam_nop_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_responds_sam_nop_500)]; markersize=10, color=:blue, marker='O', label="RBM (❌)")

	Makie.scatter!(ax, aptamer_natural_distances_rep0[rbm_seqs_rep0 ∩ findall(_inconclusive_rep0)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_inconclusive_rep0)]; markersize=15, color=(:gray, 0.5))
	Makie.scatter!(ax, aptamer_natural_distances_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_yes_rep0)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_yes_rep0)]; markersize=10, color=:blue, marker='●')
	Makie.scatter!(ax, aptamer_natural_distances_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_nop_rep0)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_nop_rep0)]; markersize=10, color=:blue, marker='O')

	Makie.ylims!(ax, 240, 360)
	Makie.axislegend(ax, position=(-0.02, -0.01), framevisible=false, nbanks=1, colgap=1, rowgap=0.1, patchlabelgap=0)

	ax = Makie.Axis(fig[1,2], width=300, height=300, xlabel="Divergence from closest natural", ylabel="RBM score", title="Distance along conserved sites", xticks=0:0.1:0.6, yticks=200:25:350, xgridvisible=false, ygridvisible=false)

	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500_conserved[rbm_seqs_500 ∩ findall(_inconclusive_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_inconclusive_500)]; markersize=15, color=(:gray, 0.5))
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500_conserved[rbm_seqs_500 ∩ findall(_responds_sam_yes_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_responds_sam_yes_500)]; markersize=10, color=:blue, marker='●')
	Makie.scatter!(ax, vec(minimum(aptamer_natural_distances_500_conserved[rbm_seqs_500 ∩ findall(_responds_sam_nop_500), :]; dims=2)) / 108, -aptamer_rbm_energies_500[rbm_seqs_500 ∩ findall(_responds_sam_nop_500)]; markersize=10, color=:blue, marker='O')

	Makie.scatter!(ax, aptamer_natural_distances_rep0_conserved[rbm_seqs_rep0 ∩ findall(_inconclusive_500)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_inconclusive_500)]; markersize=15, color=(:gray, 0.5))
	Makie.scatter!(ax, aptamer_natural_distances_rep0_conserved[rbm_seqs_rep0 ∩ findall(_responds_sam_yes_rep0)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_yes_rep0)]; markersize=10, color=:blue, marker='●')
	Makie.scatter!(ax, aptamer_natural_distances_rep0_conserved[rbm_seqs_rep0 ∩ findall(_responds_sam_nop_rep0)] / 108, -aptamer_rbm_energies_rep0[rbm_seqs_rep0 ∩ findall(_responds_sam_nop_rep0)]; markersize=10, color=:blue, marker='O')

	Makie.ylims!(ax, 240, 360)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA/cossio/SAM/2024/SamApp2025.jl/pluto/SI/Figures/Hamming distances conserved.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═b260395b-0f33-4289-83d8-8fdfa851c78c
# ╠═1dec11ae-386c-11ef-1b4c-51f221040ed7
# ╠═d43e0635-c1f2-48ff-b92f-ae7000f305a3
# ╠═1786306b-e9a5-4013-a1f6-a1133c3f99e3
# ╠═96134a6b-a910-4831-9c75-e7d052b04481
# ╠═eb41dfa5-12cc-4262-b0ce-15ed56c3a2cb
# ╠═e9a45fb3-7b6f-4957-8f8b-240c964e3fb6
# ╠═072312ba-58f2-4bea-898c-5164c976ea16
# ╠═0e9405c3-452c-4278-8b31-8aa89af09d06
# ╠═8d4bc108-2cb6-416c-9a9d-e5cdf2120206
# ╠═8798da86-73ba-4474-81dd-10296384e48e
# ╠═48a3772d-1f32-490c-9978-3099beb188dc
# ╠═9853b581-ca77-4b4d-ac8b-6ca3a960d244
# ╠═3cc2bb41-3213-42d4-995e-736eeaf01534
# ╠═26683581-756c-47c3-98bf-91a5035875d0
# ╠═4d577148-6793-493c-be1e-50d447224f47
# ╠═2f5a46ec-7ab2-444e-b050-b2c2c45a7e4d
# ╠═8d4140bf-d8cf-45fa-8ced-4bd7900a1721
# ╠═6e530ef6-e22c-4d9d-a7d3-523682526e09
# ╠═2326a0b9-2373-4629-b67e-f190f6d37102
# ╠═9483c0e6-2fa1-489e-bfdc-625b766bc2a9
# ╠═4ab7bcfc-5d88-4820-a623-e59f27550cde
# ╠═58f69627-684e-4951-a447-57240bb9e38a
# ╠═eaec262e-d68c-49f5-b14a-d4c1e6a10ee1
# ╠═3af22d1d-dc5b-45e9-953f-054320b33deb
# ╠═b5bd2ee2-289c-4ae8-8165-7ea81fb2bdfb
# ╠═fb2046e6-a7fb-4ebf-9b0b-d34865c5f9a9
# ╠═fbcd5369-f492-4857-bb30-7fa695bd06c5
# ╠═c43cf4ff-1884-427b-9ac0-d4f9efa75a8a
# ╠═33775511-874f-415a-8486-1bb317973ad0
# ╠═1b3278e0-231c-47de-ad66-068551a084c7
# ╠═609d1eea-270b-4ffd-8706-9b14849cc76a
# ╠═efa6e6c0-f018-4c4c-9bca-8267947265f1
# ╠═a4f602b8-5a0b-4147-9643-5e81aba14f48
# ╠═4f39ca7b-75ab-4168-9314-64a06d1779e7
# ╠═376848db-5682-47d2-8279-0eb8b30e1bce
# ╠═1f2365e4-f6ef-449a-b244-073b4099fb3b
# ╠═cafde511-5ef8-4684-a945-90139bb20826
# ╠═94439257-4a29-4a80-bb1c-e9afc7da3753
# ╠═612e5568-84c5-4b8b-a94d-3277c05eac0d
# ╠═5c66a219-fabf-4187-aad0-fd7ddcdd675b
# ╠═a6adc734-4075-45b1-b11b-eee8b5d91677
# ╠═de88df39-55c1-4518-8a9e-82e66b471a3d
# ╠═d28b51e3-c24d-4636-958c-cdf66fccd09b
# ╠═120eac94-3480-4ba7-ba43-15a9fb84856c
# ╠═da8d155e-bd5b-4eff-a830-378474d0d9c4
# ╠═88a27614-0463-44ee-8508-76b65e0936b9
# ╠═60ae1097-6576-40d6-a813-a2ff4f7975ec
# ╠═7ce05caa-e694-4a61-8085-7ffc3a732667
# ╠═1f090e07-cd4e-4ded-8da9-0d3be10451c2
# ╠═8a60ae44-7aad-4265-93c1-1caba00afd13
# ╠═7dcb4fff-c8ab-4f3a-b36c-19b9db528142
# ╠═8ec5db2f-da12-48f0-8e1e-7a8c97901742
# ╠═2d45af78-aa2e-4d8f-8fea-1e8e865abcb6
# ╠═bdec681b-179f-4e7e-9b5d-b2a4ca01d5af
# ╠═87f6e1da-8490-45b3-8eff-6a78bc6a7bc9
# ╠═5abad490-b3d9-4fb2-a87c-89b160c32b0b
# ╠═01af5f88-a0ed-4ee5-8c8a-b744463f4313
# ╠═b17a6262-e87e-4344-aa41-35cf78770d86
# ╠═fbac32fe-ecc7-443b-8f2a-6e0935531cba
# ╠═09d81eb5-4759-486c-8e0a-71f7f56fb7cb
# ╠═f1700a78-538c-41d1-bb60-d644e0b27cec
# ╠═b9689ab6-a819-4f22-b8bc-1add6e7366c7
# ╠═7aa22d92-c726-401f-ad37-e4efd683df6e

### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 77855f93-2b64-45bc-a307-df7f6e6187b3
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
using BioSequences: LongRNA

# ╔═╡ 55735555-c2f7-41f0-becd-cbad5717d7be
using DataFrames: DataFrame

# ╔═╡ fb8e8cbd-050d-4d0e-81ac-32528f628a0e
using Makie: @L_str

# ╔═╡ 48ab42a5-3ea8-4b69-acfb-3e2471de7144
using NaNStatistics: nanstd, nanmean, nansum, nancor

# ╔═╡ e73b3886-162a-4a77-a28b-cdf269109b98
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2abcc581-8722-4cf7-bc09-8bf98a9b8648
using Statistics: cor

# ╔═╡ 6d6e3dc8-28eb-496d-9622-8128422ed2ed
using Statistics: mean

# ╔═╡ c7378daf-3efe-4317-a95e-cd538a64e1a6
using Statistics: std

# ╔═╡ 7edc30fd-a2d6-4818-855c-c7f69c9f589b
using Statistics: middle

# ╔═╡ 4974c2e2-058d-41ca-924d-16709e4a58e6
using StatsBase: countmap

# ╔═╡ 94d99837-7415-4acb-b5d3-3b1dec5af05e
using Unitful: ustrip

# ╔═╡ f6d726bd-4493-4aee-a824-a36c72b16e95
using SamApp2025: successes_tuple_str

# ╔═╡ a06dbeed-cd34-464f-95fc-f3659f95f760
md"# Imports"

# ╔═╡ 45907f4d-29ec-4be7-97e8-bfcb4695416b
import CairoMakie

# ╔═╡ f743167e-8552-4002-b9ca-c995cf3e9829
import CSV

# ╔═╡ 6be9532a-762f-4830-8d7f-175545fe75c9
import FASTX

# ╔═╡ 02b62a0f-159c-470d-91d0-673898ca277b
import HDF5

# ╔═╡ 651a9155-c7fc-4858-ae65-b433d1d116cc
import Infernal

# ╔═╡ c6371fcc-bb84-4ae4-a2d8-972b362c5350
import KernelDensity

# ╔═╡ 297c712c-e7e6-4d2f-a9c2-59b29665e6e3
import Makie

# ╔═╡ c53d715e-40d3-44cb-b85a-9c4c61d99819
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ f513852e-ce6a-41f0-98aa-155e06244aaf
import Rfam

# ╔═╡ ca499e53-296f-4cf3-9df1-5070e22fd6f0
import SamApp2025

# ╔═╡ 15974106-a5c6-4867-acfb-cdc69f5a57be
import StatsBase

# ╔═╡ b471828b-0822-4a08-a63c-5fa01e8d90b2
import ViennaRNA

# ╔═╡ 66de62df-32fb-4dc0-b3b5-05d0e5ca718c
import PlutoUI

# ╔═╡ 02bcf2aa-204e-4385-8790-3c76432deacf
PlutoUI.TableOfContents()

# ╔═╡ b02ae72e-1892-4fdc-8c71-4e2007c65895
md"# Load data"

# ╔═╡ 2a3fff7d-a1e5-474e-8981-ac9c05494d4d
@show Rfam.get_rfam_directory() Rfam.get_rfam_version();

# ╔═╡ c182ce25-a47f-4368-9ca7-9b5b93983dcb
# load SHAPE data
shape_data_rep123 = SamApp2025.load_shapemapper_data_pierre_demux_20230929_repl123();

# ╔═╡ 46c07093-3c11-45b8-b120-c9ed429bc458
# load SHAPE data
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 3105982c-f2a4-4f0a-b09c-29a7cf76f200
# split rep0 from rep4+5
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 3c96e229-e909-4de6-849e-753109b229d5
shape_data_rep123.conditions

# ╔═╡ f2a92753-4792-41d5-a788-972c8f03e7c5
conds_mg_rep123 = findall(==("SAMAP_1M7_noSAM_5Mg_T30C_rep6"),  shape_data_rep123.conditions)

# ╔═╡ 775e36c0-7c9b-47e7-adef-ac686ce0dc88
conds_sam_rep123 = findall(∈(["SAMAP_1M7_1SAM_5Mg_T30C_rep6", "SAMAP_1M7_0-1SAM_5Mg_T30C_rep6"]), shape_data_rep123.conditions)

# ╔═╡ e1580d5b-8d48-43bb-8da2-b723f04582f5
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ a4f61bcd-441e-4f2b-b6cc-0f9e997be95c
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 8759ae29-c95f-4bf4-b1a0-fb58689550c5
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ ba2ea96a-aa4c-467d-82b5-688b46cc67bd
(; bps, nps, pks) = SamApp2025.RF00162_sites_paired()

# ╔═╡ f74e788f-47a3-46d6-80da-8e3e8b50a029
full_seqs_rep123 = findall(shape_data_rep123.aptamer_origin .== "RF00162_full30")

# ╔═╡ 802a5655-9b21-448b-933d-0471350f1058
seed_seqs_rep123 = findall(shape_data_rep123.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 9504ca2a-adad-4a04-a943-9b4cee91f834
nat_seqs_rep123 = full_seqs_rep123 ∪ seed_seqs_rep123;

# ╔═╡ e77f6789-83be-4d6e-8b58-30e2ac154777
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 0e13e6a9-4fad-495c-b6c9-c6412967da55
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 0c468817-a3d5-4711-941a-c5616ad0c662
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ e5230e7c-4c61-457f-b57c-07a5272074b5
bps_reactivities_rep123 = shape_data_rep123.shape_reactivities[bps, nat_seqs_rep123, conds_sam_rep123];

# ╔═╡ fcd6da1f-2afb-4b47-8a39-4aacb2b48fb2
nps_reactivities_rep123 = shape_data_rep123.shape_reactivities[nps, nat_seqs_rep123, conds_sam_rep123];

# ╔═╡ cb886986-01e2-494b-8d5d-0ceacaa6d886
all_reactivities_rep123 = shape_data_rep123.shape_reactivities[:, nat_seqs_rep123, conds_sam_rep123];

# ╔═╡ f1546c98-b0ad-45f9-b3a6-4eebc5a4160a
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ b6eb5309-f55f-426a-ab1d-45b3717d9782
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 9dba8ed2-2be4-4b24-87ec-4dd94f3aa286
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 92975e4a-13a1-42b2-8219-a1266b7189a1
shape_stats_rep123 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep123,
    paired_reactivities = bps_reactivities_rep123,
    unpaired_reactivities = nps_reactivities_rep123,
    all_reactivities = all_reactivities_rep123,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 427e0bce-110b-4b8c-af80-d2f0f9f38503
shape_stats_rep0 = SamApp2025.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ ac0f6b9c-abdc-4aa3-bf00-29f63330d219
_thresh = log(5)

# ╔═╡ 3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
_sites = SamApp2025.hallmark_sites_20230507;

# ╔═╡ acf8bb76-bf87-4cf8-b280-e0e13c65fdff
x_mg_rep123 = nansum(shape_stats_rep123.shape_log_odds[_sites, :,  conds_mg_rep123]; dim=(1,3))

# ╔═╡ 1a364c8f-49e3-4c33-b892-8f35106c4afa
x_sam_rep123 = nansum(shape_stats_rep123.shape_log_odds[_sites, :, conds_sam_rep123]; dim=(1,3))

# ╔═╡ d71e1b09-2f2b-40ab-a52d-8bc43bfd6403
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))

# ╔═╡ 70dcf640-9822-4175-8275-03000b86862c
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

# ╔═╡ 4b72ca7b-1d8a-4509-aef8-026b2d9b644a
_responds_sam_yes_rep123 = (x_mg_rep123 .< -_thresh) .& (x_sam_rep123 .> +_thresh);

# ╔═╡ 11b9e4bf-1531-459c-9d36-6d73526a96fb
_responds_sam_nop_rep123 = (x_mg_rep123 .> +_thresh) .| (x_sam_rep123 .< -_thresh);

# ╔═╡ 4825347b-2117-4b48-9f72-3bc24be76d07
_inconclusive_rep123 = ((!).(_responds_sam_yes_rep123)) .& ((!).(_responds_sam_nop_rep123));

# ╔═╡ 8f4dbcc1-3f77-43f0-b4bb-bb027bf20569
_conclusive_rep123 = _responds_sam_yes_rep123 .| _responds_sam_nop_rep123;

# ╔═╡ 23069178-e9c4-448f-95b0-45ed23411926
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 8b90364d-67cd-4ca6-bba9-062d5c1e409e
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ 64cb1776-0d3f-49ca-9b3f-f27b3bcb9acb
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ f21cff04-7336-4ef9-9358-5a5c1c8a2c76
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

# ╔═╡ db615f71-8858-47da-be2f-8196f7070e6d
sum(_responds_sam_yes_rep123), sum(_responds_sam_nop_rep123), sum(_inconclusive_rep123)

# ╔═╡ 41caae29-9181-4cf8-a814-64868a9b8e90
aptamer_rbm_energies = [
    ismissing(seq) ? missing :
    free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(LongRNA{4}(seq)))
    for seq in shape_data_rep123.aligned_sequences
];

# ╔═╡ 3d0cd846-07a7-4910-a2a2-2ef660355ad5
wuss = SamApp2025.rfam_ss("RF00162"; inserts=false)

# ╔═╡ 74fbb354-2657-475e-8d1e-42d39e255396
ss = SamApp2025.clean_wuss(wuss)

# ╔═╡ 9f0baca9-46e2-4b20-9c0c-0dbbac36bbf0
p1_pos = SamApp2025.RF00162_sites_annotated_secondary_structure().p1;

# ╔═╡ ffcd720b-134f-4b83-adf7-ab39ab6dfba1
p2_pos = SamApp2025.RF00162_sites_annotated_secondary_structure().p2;

# ╔═╡ de174fc5-5bc0-4b5c-980d-88168f5355d1
p3_pos = SamApp2025.RF00162_sites_annotated_secondary_structure().p3;

# ╔═╡ bea93ba6-a1e5-444a-b419-f0ec290cb61f
p4_pos = SamApp2025.RF00162_sites_annotated_secondary_structure().p4;

# ╔═╡ 2639031c-df64-47c4-9d75-4f2ad73291dd
pk_pos = SamApp2025.RF00162_sites_annotated_secondary_structure().pk;

# ╔═╡ 15ff286a-a8e4-4c8b-9c73-04450ad6c500
ss_without_P1 = join([i ∈ p1_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ b798af09-363d-4800-b37e-a8fe8ce36398
ss_without_P2 = join([i ∈ p2_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 1e382321-089d-407f-9d2b-e07dcdb00135
ss_without_P3 = join([i ∈ p3_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 445ab169-1024-4238-b3bb-35b9cbdc99e7
ss_without_P4 = join([i ∈ p4_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 6212efdb-8f9a-4fc8-824c-0623ac71d2ce
ss_pk_only = replace(wuss, r"\(|\)|\[|\]|\{|\}|\<|\>|\-|\_|\," => '.', 'A' => '(', 'a' => ')')

# ╔═╡ 34cc67fd-f567-425b-bad0-117b00c85a1d
md"# Counts"

# ╔═╡ dd25496f-e430-42df-a0b4-143bcb26f9f3
begin
	println("Replicate (naturals)")

	println("Nat. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[nat_seqs_rep123]),
	        sum(_responds_sam_nop_rep123[nat_seqs_rep123])
	    )
	)
	println("Nat seed. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[seed_seqs_rep123]),
	        sum(_responds_sam_nop_rep123[seed_seqs_rep123])
	    )
	)
	println("Nat full. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[full_seqs_rep123]),
	        sum(_responds_sam_nop_rep123[full_seqs_rep123])
	    )
	)
	println("Nat. (E<-300): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ nat_seqs_rep123]),
	        sum(_responds_sam_nop_rep123[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ nat_seqs_rep123])
	    )
	)
	println("Nat. (E<-310): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ nat_seqs_rep123]),
	        sum(_responds_sam_nop_rep123[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ nat_seqs_rep123])
	    )
	)

	println("All. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123),
	        sum(_responds_sam_nop_rep123)
	    )
	)
	println("All. (E<-300): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[findall(replace(aptamer_rbm_energies .< -300, missing => false))]),
	        sum(_responds_sam_nop_rep123[findall(replace(aptamer_rbm_energies .< -300, missing => false))])
	    )
	)
	println("All. (E<-310): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep123[findall(replace(aptamer_rbm_energies .< -310, missing => false))]),
	        sum(_responds_sam_nop_rep123[findall(replace(aptamer_rbm_energies .< -310, missing => false))])
	    )
	)
end

# ╔═╡ eff6428b-4022-4d8a-8cd0-19674d8dadcb
_responds_sam_yes_rep123[nat_seqs_rep123]

# ╔═╡ 30fb8f07-b981-491e-a0ec-fa3b47ef1fc9
mean(_responds_sam_yes_rep0[nat_seqs_rep0] .== _responds_sam_yes_rep123[nat_seqs_rep123])

# ╔═╡ c0e7239d-21c0-4cc4-9cce-c9870b64574d
length(findall(_responds_sam_yes_rep0[nat_seqs_rep0]) ∩ findall(_responds_sam_yes_rep123[nat_seqs_rep123])) / sum(_responds_sam_yes_rep0[nat_seqs_rep0])

# ╔═╡ 4e02396a-8314-41d9-b7de-e893ed196fa1
length(findall(_responds_sam_yes_rep0[nat_seqs_rep0]) ∩ findall(_responds_sam_yes_rep123[nat_seqs_rep123])) / sum(_responds_sam_yes_rep123[nat_seqs_rep123])

# ╔═╡ 1c958bae-4df5-4df3-ba59-a9bfb0bc2e49
length(findall(_responds_sam_yes_rep0[nat_seqs_rep0]) ∩ findall(_responds_sam_yes_rep123[nat_seqs_rep123]))

# ╔═╡ b3c1f337-669f-4f3f-8787-ecdc925c3b2e
sum(_responds_sam_yes_rep123[nat_seqs_rep123])

# ╔═╡ 55b64786-b444-45e6-ba46-fb13a20216f8
length(findall(_responds_sam_nop_rep0[nat_seqs_rep0]) ∩ findall(_responds_sam_nop_rep123[nat_seqs_rep123])) / sum(_responds_sam_nop_rep123[nat_seqs_rep123])

# ╔═╡ e2e43989-0eea-4faf-8b6c-a779f71f9dd3
sum(_responds_sam_yes_rep0[nat_seqs_rep0])

# ╔═╡ 39db31d8-e978-4c10-91da-cd997b3f07e0
length(_responds_sam_yes_rep123)

# ╔═╡ 29989b1f-6d1a-4b3f-8b3f-88394e26a37e
md"# Replicate correlations"

# ╔═╡ 54d12444-8644-4667-b7d3-0d7a8ea6a321
let fig = Makie.Figure()
	_sz = 300
	_err_thresh = 0.5

	# naturals only

	_x = vec(shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_noSAM_5Mg_T30C_rep0"))])
	_y = vec(shape_data_rep123.shape_reactivities[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_1SAM_5Mg_T30C_rep6"))])
	_x_err = vec(shape_data_rep0.shape_reactivities_err[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_noSAM_5Mg_T30C_rep0"))])
	_y_err = vec(shape_data_rep123.shape_reactivities_err[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_noSAM_5Mg_T30C_rep6"))])
	_idx_finite = findall(isfinite, _x) ∩ findall(isfinite, _y)
	_idx_filter = findall(_x_err .< _err_thresh * sqrt.(abs.(_x))) ∩ findall(_y_err .< _err_thresh * sqrt.(abs.(_y)))

	ax = Makie.Axis(fig[1,1], width=_sz, height=_sz, xlabel="Reactivity (Repl.0)", ylabel="Reactivity (Repl.1)", title="no SAM + 5Mg", xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, _x[_idx_finite], _y[_idx_finite], markersize=10, color=(:gray, 0.5), label="all (cor. = $(round(cor(_x[_idx_finite], _y[_idx_finite]), digits=2)))")
	Makie.scatter!(ax, _x[_idx_filter], _y[_idx_filter], markersize=5, color=:blue, label="filtered (cor. = $(round(cor(_x[_idx_filter], _y[_idx_filter]), digits=2)))")
	Makie.axislegend(ax; framevisible=false, position=:lt)


	_x = vec(shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_0-1SAM_5Mg_T30C_rep0"))])
	_y = vec(shape_data_rep123.shape_reactivities[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_0-1SAM_5Mg_T30C_rep6"))])
	_x_err = vec(shape_data_rep0.shape_reactivities_err[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_0-1SAM_5Mg_T30C_rep0"))])
	_y_err = vec(shape_data_rep123.shape_reactivities_err[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_0-1SAM_5Mg_T30C_rep6"))])
	_idx = findall(_x_err .< _err_thresh * sqrt.(abs.(_x))) ∩ findall(_y_err .< _err_thresh * sqrt.(abs.(_y)))

	ax = Makie.Axis(fig[1,2], width=_sz, height=_sz, xlabel="Reactivity (Repl.0)", ylabel="Reactivity (Repl.1)", title="0.1SAM + 5Mg", xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, _x, _y, markersize=10, color=(:gray, 0.5), label="all (Pearson = $(round(nancor(_x, _y), digits=2)))")
	Makie.scatter!(ax, _x[_idx], _y[_idx], markersize=5, color=:blue, label="filtered (Pearson = $(round(nancor(_x[_idx], _y[_idx]), digits=2)))")
	Makie.axislegend(ax; framevisible=false, position=:lt)


	_x = vec(shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_1SAM_5Mg_T30C_rep0"))])
	_y = vec(shape_data_rep123.shape_reactivities[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_1SAM_5Mg_T30C_rep6"))])
	_x_err = vec(shape_data_rep0.shape_reactivities_err[:, nat_seqs_rep0, only(findall(shape_data_rep0.conditions .== "SAMAP_1M7_1SAM_5Mg_T30C_rep0"))])
	_y_err = vec(shape_data_rep123.shape_reactivities_err[:, nat_seqs_rep123, only(findall(shape_data_rep123.conditions .== "SAMAP_1M7_1SAM_5Mg_T30C_rep6"))])
	_idx = findall(_x_err .< _err_thresh * sqrt.(abs.(_x))) ∩ findall(_y_err .< _err_thresh * sqrt.(abs.(_y)))

	ax = Makie.Axis(fig[1,3], width=_sz, height=_sz, xlabel="Reactivity (Repl.0)", ylabel="Reactivity (Repl.1)", title="1SAM + 5Mg", xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, _x, _y, markersize=10, color=(:gray, 0.5), label="all (Pearson = $(round(nancor(_x, _y), digits=2)))")
	Makie.scatter!(ax, _x[_idx], _y[_idx], markersize=5, color=:blue, label="filtered (Pearson = $(round(nancor(_x[_idx], _y[_idx]), digits=2)))")
	Makie.axislegend(ax; framevisible=false, position=:lt)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Repl0_vs_Repl123_reactivities.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═a06dbeed-cd34-464f-95fc-f3659f95f760
# ╠═77855f93-2b64-45bc-a307-df7f6e6187b3
# ╠═45907f4d-29ec-4be7-97e8-bfcb4695416b
# ╠═f743167e-8552-4002-b9ca-c995cf3e9829
# ╠═6be9532a-762f-4830-8d7f-175545fe75c9
# ╠═02b62a0f-159c-470d-91d0-673898ca277b
# ╠═651a9155-c7fc-4858-ae65-b433d1d116cc
# ╠═c6371fcc-bb84-4ae4-a2d8-972b362c5350
# ╠═297c712c-e7e6-4d2f-a9c2-59b29665e6e3
# ╠═c53d715e-40d3-44cb-b85a-9c4c61d99819
# ╠═f513852e-ce6a-41f0-98aa-155e06244aaf
# ╠═ca499e53-296f-4cf3-9df1-5070e22fd6f0
# ╠═15974106-a5c6-4867-acfb-cdc69f5a57be
# ╠═b471828b-0822-4a08-a63c-5fa01e8d90b2
# ╠═66de62df-32fb-4dc0-b3b5-05d0e5ca718c
# ╠═c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
# ╠═55735555-c2f7-41f0-becd-cbad5717d7be
# ╠═fb8e8cbd-050d-4d0e-81ac-32528f628a0e
# ╠═48ab42a5-3ea8-4b69-acfb-3e2471de7144
# ╠═e73b3886-162a-4a77-a28b-cdf269109b98
# ╠═2abcc581-8722-4cf7-bc09-8bf98a9b8648
# ╠═6d6e3dc8-28eb-496d-9622-8128422ed2ed
# ╠═c7378daf-3efe-4317-a95e-cd538a64e1a6
# ╠═7edc30fd-a2d6-4818-855c-c7f69c9f589b
# ╠═4974c2e2-058d-41ca-924d-16709e4a58e6
# ╠═94d99837-7415-4acb-b5d3-3b1dec5af05e
# ╠═f6d726bd-4493-4aee-a824-a36c72b16e95
# ╠═02bcf2aa-204e-4385-8790-3c76432deacf
# ╠═b02ae72e-1892-4fdc-8c71-4e2007c65895
# ╠═2a3fff7d-a1e5-474e-8981-ac9c05494d4d
# ╠═c182ce25-a47f-4368-9ca7-9b5b93983dcb
# ╠═46c07093-3c11-45b8-b120-c9ed429bc458
# ╠═3105982c-f2a4-4f0a-b09c-29a7cf76f200
# ╠═3c96e229-e909-4de6-849e-753109b229d5
# ╠═f2a92753-4792-41d5-a788-972c8f03e7c5
# ╠═775e36c0-7c9b-47e7-adef-ac686ce0dc88
# ╠═e1580d5b-8d48-43bb-8da2-b723f04582f5
# ╠═a4f61bcd-441e-4f2b-b6cc-0f9e997be95c
# ╠═8759ae29-c95f-4bf4-b1a0-fb58689550c5
# ╠═ba2ea96a-aa4c-467d-82b5-688b46cc67bd
# ╠═f74e788f-47a3-46d6-80da-8e3e8b50a029
# ╠═802a5655-9b21-448b-933d-0471350f1058
# ╠═9504ca2a-adad-4a04-a943-9b4cee91f834
# ╠═e77f6789-83be-4d6e-8b58-30e2ac154777
# ╠═0e13e6a9-4fad-495c-b6c9-c6412967da55
# ╠═0c468817-a3d5-4711-941a-c5616ad0c662
# ╠═e5230e7c-4c61-457f-b57c-07a5272074b5
# ╠═fcd6da1f-2afb-4b47-8a39-4aacb2b48fb2
# ╠═cb886986-01e2-494b-8d5d-0ceacaa6d886
# ╠═f1546c98-b0ad-45f9-b3a6-4eebc5a4160a
# ╠═b6eb5309-f55f-426a-ab1d-45b3717d9782
# ╠═9dba8ed2-2be4-4b24-87ec-4dd94f3aa286
# ╠═92975e4a-13a1-42b2-8219-a1266b7189a1
# ╠═427e0bce-110b-4b8c-af80-d2f0f9f38503
# ╠═ac0f6b9c-abdc-4aa3-bf00-29f63330d219
# ╠═3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
# ╠═acf8bb76-bf87-4cf8-b280-e0e13c65fdff
# ╠═1a364c8f-49e3-4c33-b892-8f35106c4afa
# ╠═d71e1b09-2f2b-40ab-a52d-8bc43bfd6403
# ╠═70dcf640-9822-4175-8275-03000b86862c
# ╠═4b72ca7b-1d8a-4509-aef8-026b2d9b644a
# ╠═11b9e4bf-1531-459c-9d36-6d73526a96fb
# ╠═4825347b-2117-4b48-9f72-3bc24be76d07
# ╠═8f4dbcc1-3f77-43f0-b4bb-bb027bf20569
# ╠═23069178-e9c4-448f-95b0-45ed23411926
# ╠═8b90364d-67cd-4ca6-bba9-062d5c1e409e
# ╠═64cb1776-0d3f-49ca-9b3f-f27b3bcb9acb
# ╠═f21cff04-7336-4ef9-9358-5a5c1c8a2c76
# ╠═db615f71-8858-47da-be2f-8196f7070e6d
# ╠═41caae29-9181-4cf8-a814-64868a9b8e90
# ╠═3d0cd846-07a7-4910-a2a2-2ef660355ad5
# ╠═74fbb354-2657-475e-8d1e-42d39e255396
# ╠═9f0baca9-46e2-4b20-9c0c-0dbbac36bbf0
# ╠═ffcd720b-134f-4b83-adf7-ab39ab6dfba1
# ╠═de174fc5-5bc0-4b5c-980d-88168f5355d1
# ╠═bea93ba6-a1e5-444a-b419-f0ec290cb61f
# ╠═2639031c-df64-47c4-9d75-4f2ad73291dd
# ╠═15ff286a-a8e4-4c8b-9c73-04450ad6c500
# ╠═b798af09-363d-4800-b37e-a8fe8ce36398
# ╠═1e382321-089d-407f-9d2b-e07dcdb00135
# ╠═445ab169-1024-4238-b3bb-35b9cbdc99e7
# ╠═6212efdb-8f9a-4fc8-824c-0623ac71d2ce
# ╠═34cc67fd-f567-425b-bad0-117b00c85a1d
# ╠═dd25496f-e430-42df-a0b4-143bcb26f9f3
# ╠═eff6428b-4022-4d8a-8cd0-19674d8dadcb
# ╠═30fb8f07-b981-491e-a0ec-fa3b47ef1fc9
# ╠═c0e7239d-21c0-4cc4-9cce-c9870b64574d
# ╠═4e02396a-8314-41d9-b7de-e893ed196fa1
# ╠═1c958bae-4df5-4df3-ba59-a9bfb0bc2e49
# ╠═b3c1f337-669f-4f3f-8787-ecdc925c3b2e
# ╠═55b64786-b444-45e6-ba46-fb13a20216f8
# ╠═e2e43989-0eea-4faf-8b6c-a779f71f9dd3
# ╠═39db31d8-e978-4c10-91da-cd997b3f07e0
# ╠═29989b1f-6d1a-4b3f-8b3f-88394e26a37e
# ╠═54d12444-8644-4667-b7d3-0d7a8ea6a321

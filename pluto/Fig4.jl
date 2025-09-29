### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 952a79aa-12fa-11ef-2326-d5fe00ff3dfe
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ adbed31f-5777-4b6f-acf9-f37a472ef5d8
using BioSequences: LongRNA

# ╔═╡ de469256-6471-4fa9-854b-7b7783295aef
using DataFrames: DataFrame

# ╔═╡ a9fb22f9-4c47-4430-a0aa-cc4075b86f1f
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ ebabc892-15f1-4b91-be5d-b85676d94fa7
using LinearAlgebra: Diagonal, eigen

# ╔═╡ 2a184e1f-bcba-48a3-a54b-793a25d9501e
using Makie: @L_str

# ╔═╡ 05b55f01-751c-4910-964f-aadb33f844c0
using NaNStatistics: nansum

# ╔═╡ a72c9391-c8f5-4bf5-94bd-3e9bbc326fa2
using Random: bitrand

# ╔═╡ bf22e1b2-f8bd-4990-938d-4db503660c3d
using Statistics: cor, mean

# ╔═╡ ac2f96ec-b31f-4f63-8fd3-4390d687755b
using StatsBase: countmap

# ╔═╡ 5d5eb95f-6801-47c3-a82e-afab35559c0e
md"# Imports"

# ╔═╡ f01ccb33-26f8-4fb9-94a2-88be47db3762
import Makie, CairoMakie, PlutoUI

# ╔═╡ 897ffcd7-e59a-4121-975d-29bc9a27988b
import CSV, HDF5

# ╔═╡ 49b0a181-8c3b-41dc-b3e0-9e55efa87bd9
import FASTX, Infernal, Rfam

# ╔═╡ 4a2a79bc-6fb1-45d3-ad4a-abfc889f2737
import SamApp2025

# ╔═╡ 5ca1d179-8dda-4e4b-9bab-b1cf4f6bafbe
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 98c97d8b-9849-48dc-ae6c-c0924d14e560
import StatsBase, KernelDensity

# ╔═╡ aa4f2963-ec25-4469-adea-f7f16c69db26
PlutoUI.TableOfContents()

# ╔═╡ aa8a002e-693a-402b-aaf1-576bf59e3a4e
md"# Plot"

# ╔═╡ 3a00d0fd-17e1-4518-8fca-ac1ccd0f7d28
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162")

# ╔═╡ 660fa41b-c6c2-4b41-a698-82530c8dc4ef
RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")

# ╔═╡ 8399bdac-bd99-4a20-be2c-2055db32a5a6
RF00162_seed_match_cols = findall(≠('.'), SamApp2025.stockholm_ss(RF00162_seed_stk.out));

# ╔═╡ fd9539ff-e725-4b18-94fd-e6a75fc2ec73
RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM") # WARNING: this has inserts marked as '-'

# ╔═╡ 6538d3e8-54d8-484d-8baa-3d3cc44ed72c
RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))

# ╔═╡ 71b73684-6053-4efc-92b7-391267f26ad9
RF00162_seed_seqs_noinserts = LongRNA{4}.([FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records]);

# ╔═╡ e4309e31-4086-4b3f-baec-8798ab975c44
# trimmed (no inserts) aligned fasta

# ╔═╡ 3a5e2efa-a739-42c0-a76c-50c6c07487e9
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");
# these are already aligned and without inserts

# ╔═╡ 0d4f2776-4666-48dd-80df-990c479569df
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ 488c249d-9810-4697-a6ca-7ab326b53762
# emit sequences from Rfam CM model
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Rfam_cm.out; N=5000, aligned=true, outformat="AFA");

# ╔═╡ 0eeacfc3-2240-46e2-b851-9d574771d1d7
begin
	Rfam_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)));
	Rfam_cm_emitted_sequences = [filter(!=('.'), filter(!islowercase, seq)) for seq in Rfam_cm_emitted_sequences];
	Rfam_cm_emitted_sequences = LongRNA{4}.(Rfam_cm_emitted_sequences);
end

# ╔═╡ 896d95be-2a0c-4c4b-9c73-d6e2d1452306
begin
	# aligned hits, used to train a new noiseless CM model (in Stockholm format, without inserts!)
	RF00162_hits_stk = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true);
	# fit new CM model using full alignment (without inserts), and without entropic noise
	Refined_cm = Infernal.cmbuild(RF00162_hits_stk.out; enone=true);

	# emit sequences from Refined CM model
	Refined_cm_emitted_sequences_afa = Infernal.cmemit(Refined_cm.cmout; N=5000, aligned=true, outformat="AFA");
	Refined_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Refined_cm_emitted_sequences_afa.out)));
	# remove inserts
	Refined_cm_emitted_sequences = [filter(!=('.'), filter(!islowercase, seq)) for seq in Refined_cm_emitted_sequences];
	@assert only(unique(length.(Refined_cm_emitted_sequences))) == 108
	Refined_cm_emitted_sequences = LongRNA{4}.(Refined_cm_emitted_sequences);
end

# ╔═╡ 0bd3ffae-c009-4740-aa57-4dd5126eb848
# use saved RBM samples
sampled_v = SamApp2025.rbm2022samples(); # faster

# ╔═╡ 7aba8c42-93a8-4528-90cd-b9a308960e80
begin
	# Infernal scores of hits, using Rfam CM model
	RF00162_hits_Rfam_cm_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(
			Rfam_cm.out,
			Infernal.esl_reformat("FASTA", RF00162_hits_afa.out; informat="AFA").out;
			glob=true, informat="FASTA"
		).sfile
	).bit_sc;

	# Infernal scores of hits, using Refined CM model
	RF00162_hits_Refined_cm_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(
			Refined_cm.cmout,
			Infernal.esl_reformat("FASTA", RF00162_hits_afa.out; informat="AFA").out;
			glob=true, informat="FASTA"
		).sfile
	).bit_sc;
end

# ╔═╡ 16508ca2-1c56-4136-842d-397ab4760586
begin
	# Infernal scores of Refined CM samples
	_tmpfasta = tempname()
	FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
	    for (n, seq) in enumerate(Refined_cm_emitted_sequences)
	        ismissing(seq) && continue
	        write(writer, FASTX.FASTA.Record(string(n), filter(!=('-'), string(seq))))
	    end
	end
	# Infernal scores
	Refined_cm_emitted_sequences_infernal_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(Refined_cm.cmout, _tmpfasta; glob=true, informat="FASTA").sfile
	).bit_sc;


	# Infernal scores of Rfam CM samples
	_tmpfasta = tempname()
	FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
	    for (n, seq) in enumerate(Rfam_cm_emitted_sequences)
	        ismissing(seq) && continue
	        write(writer, FASTX.FASTA.Record(string(n), filter(!=('-'), string(seq))))
	    end
	end
	# Infernal scores
	Rfam_cm_emitted_sequences_infernal_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(Rfam_cm.out, _tmpfasta; glob=true, informat="FASTA").sfile
	).bit_sc;

	# Infernal scores of RBM samples
	_tmpfasta = tempname()
	FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
	    for (n, seq) in enumerate(SamApp2025.rnaseq(sampled_v))
	        @assert !ismissing(seq)
	        write(writer, FASTX.FASTA.Record(string(n), filter(!=('-'), string(seq))))
	    end
	end

	# Rfam CM
	RBM_samples_Rfam_CM_infernal_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(Rfam_cm.out, _tmpfasta; glob=true, informat="FASTA").sfile
	).bit_sc;

	# Refined CM
	RBM_samples_Refined_CM_infernal_scores = Infernal.cmalign_parse_sfile(
		Infernal.cmalign(Refined_cm.cmout, _tmpfasta; glob=true, informat="FASTA").sfile
	).bit_sc;

end

# ╔═╡ 54d09468-5516-4568-9582-8172cc40c4f5
# sites that have some non-zero fluctuations
# We need to separate frozen sites below because otherwise cor and eigen give NaN, infinities, and fail
_variable_sites_flag = vec(all(0 .< mean(SamApp2025.onehot(RF00162_hits_sequences); dims=3) .< 1; dims=1));

# ╔═╡ 1b1989e8-36d2-474a-aa6e-be02bebb8e00
_variable_sites = findall(_variable_sites_flag);

# ╔═╡ 219bd907-7417-4718-aaf8-4555fb8f2e53
RF00162_hits_var_sites_only = SamApp2025.onehot(RF00162_hits_sequences)[:, _variable_sites, :];

# ╔═╡ cf5a9319-b99d-4c0a-a75c-f633fad52c63
RF00162_hits_cor = cor(reshape(RF00162_hits_var_sites_only, :, size(RF00162_hits_var_sites_only, 3)); dims=2);

# ╔═╡ 9f9e4969-bfbb-4f2e-971b-1265bf461a16
RF00162_hits_eig = eigen(RF00162_hits_cor);

# ╔═╡ 92fec252-d3a6-429e-ba3d-01139abc7bc7
# remap the variable sites eigenvectors back to the original consensus sequence numbering
begin
	RF00162_hits_eigvec = zeros(5, 108, size(RF00162_hits_eig.vectors, 1))
	for n in 1:size(RF00162_hits_eig.vectors, 1)
	    vec(view(RF00162_hits_eigvec, :, _variable_sites, n)) .= RF00162_hits_eig.vectors[:, n]
	end
end

# ╔═╡ 66d621fa-231f-415b-aa0b-ea6dc88cbaf1
__proj_hits = reshape(SamApp2025.onehot(RF00162_hits_sequences), 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ 68768649-9e5b-48e7-a33a-c50e20389a94
__proj_rbm = reshape(sampled_v, 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ 51789884-f149-4562-9656-06b1e210fb43
__proj_refined_cm = reshape(SamApp2025.onehot(Refined_cm_emitted_sequences), 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ b0df20e8-cd4d-4800-a348-d7bbd3e3cd65
__proj_rfam_cm = reshape(SamApp2025.onehot(Rfam_cm_emitted_sequences), 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ 93a22e9b-5f2f-4d1f-b693-e4489c6befa1
# load SHAPE data
shape_data_045 = SamApp2025.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 703fde1e-33b2-4f66-b0b7-9e883f6470d4
shape_data_rep0 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 4e1d8cc4-e169-484a-9c08-ed56db0350f5
shape_data_rep45 = SamApp2025.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ ac689699-10a5-4989-bd27-46a0fae25be2
_idx_not_missing_seqs = findall(!ismissing, shape_data_rep0.aligned_sequences);

# ╔═╡ 373b6eb8-7e76-430a-8005-6eb688ba660e
shape_sequences_onehot = SamApp2025.onehot(LongRNA{4}.(shape_data_rep0.aligned_sequences[_idx_not_missing_seqs]));

# ╔═╡ d056e7a9-3cc6-4179-9ce8-f8c567b1a56f
__proj_probed = reshape(shape_sequences_onehot, 5*108, :)' * reshape(RF00162_hits_eigvec, 5*108, :);

# ╔═╡ fd502d5b-c54d-4f02-88bb-bfacdd6e350c
_probed_origin = shape_data_rep0.aptamer_origin[_idx_not_missing_seqs];

# ╔═╡ d89ee7c9-f646-4832-bf32-427a9aa8a640
begin
	hits_tax_df = SamApp2025.rf00162_hits_taxonomy();
	hits_tax_df.taxonomy_split = [ismissing(tax) ? missing : split(tax, "; ") for tax in hits_tax_df.taxonomy]
	hits_tax_df.taxa_1 = [ismissing(tax) ? missing : length(tax) ≥ 1 ? tax[1] : missing for tax in hits_tax_df.taxonomy_split];
	hits_tax_df.taxa_2 = [ismissing(tax) ? missing : length(tax) ≥ 2 ? tax[2] : missing for tax in hits_tax_df.taxonomy_split];
	hits_tax_df.taxa_3 = [ismissing(tax) ? missing : length(tax) ≥ 3 ? tax[3] : missing for tax in hits_tax_df.taxonomy_split];
end

# ╔═╡ 38c55494-e38a-4358-bc30-ba0f2880eb10
hits_tax_cnt = countmap(split(join(filter(!ismissing, filter(!ismissing, hits_tax_df.taxonomy)), "; "), "; "));

# ╔═╡ ea3ad48b-b0ac-4054-9548-fa9ed1842a16
let fig = Makie.Figure()

	ax = Makie.Axis(fig[1,1][1,1], xlabel="rCM score", ylabel="RBM score", width=400, height=400, xticks=20:40:130, xgridvisible=false, ygridvisible=false)
	Makie.hlines!(ax, 300, color=:orange, linestyle=:dash, linewidth=2)
	#Makie.scatter!(ax, RF00162_hits_Refined_cm_scores, -RBMs.free_energy(SamApp.rbm2022(), SamApp.onehot(RF00162_hits_sequences)), label="Natural", color=(:gray, 0.5), markersize=10)
	Makie.scatter!(
		ax, RF00162_hits_Rfam_cm_scores, -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(RF00162_hits_sequences)), label="MSA", color=(:gray, 0.5), markersize=10
	)
	Makie.scatter!(ax,
	    #Refined_cm_emitted_sequences_infernal_scores[1:2000],
	    Rfam_cm_emitted_sequences_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(Refined_cm_emitted_sequences))[1:2000],
	    label="rCM", color=:red, markersize=5
	)
	Makie.scatter!(ax,
	    #RBM_samples_Refined_CM_infernal_scores[1:2000],
	    RBM_samples_Rfam_CM_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), sampled_v)[1:2000],
	    label="RBM", color=:blue, markersize=5
	)
	Makie.xlims!(ax, -7, 101)
	Makie.ylims!(ax, 220, 365)
	Makie.axislegend(ax, position=:lt, framevisible=false)


	_colors = [:purple, :cyan, :lime, :teal, :orange]
	_c = 0

	# Natural
	ax = Makie.Axis(fig[1,1][1,2], xlabel="PC1", ylabel="PC2", width=400, height=400, xgridvisible=false, ygridvisible=false) #title="Natural sequences")
	Makie.scatter!(ax, __proj_hits[:, end], __proj_hits[:, end - 1], markersize=10, label="MSA", color=(:gray, 0.5))
	for t = unique(hits_tax_df.taxa_2)
	    ismissing(t) && continue
	    hits_tax_cnt[t] > 100 || continue
	    _c += 1
	    Makie.scatter!(ax,
	        __proj_hits[replace(hits_tax_df.taxa_2 .== t, missing => false), end],
	        __proj_hits[replace(hits_tax_df.taxa_2 .== t, missing => false), end - 1],
	        markersize=5, label=t, color=_colors[_c])
	end
	Makie.axislegend(ax, position=(-0.05, -0.03), framevisible=false)

	ax = Makie.Axis(fig[2,1][1,1], xlabel="PC1", ylabel="PC2", width=250, height=250, xgridvisible=false, ygridvisible=false, title="RBM generated sequences")
	Makie.scatter!(ax, __proj_hits[:, end], __proj_hits[:, end - 1], markersize=10, label="MSA", color=color=(:gray, 0.5))
	Makie.scatter!(ax, __proj_rbm[:, end], __proj_rbm[:, end - 1], markersize=4, label="RBM", color=:blue)
	Makie.axislegend(ax, position=:lb, framevisible=false)

	ax = Makie.Axis(fig[2,1][1,2], xlabel="PC1", ylabel="PC2", width=250, height=250, xgridvisible=false, ygridvisible=false, title="CM generated sequences")
	Makie.scatter!(ax, __proj_hits[:, end], __proj_hits[:, end - 1], markersize=10, label="MSA", color=color=(:gray, 0.5))
	Makie.scatter!(ax, __proj_rfam_cm[:, end], __proj_rfam_cm[:, end - 1], markersize=4, label="rCM", color=:red)
	Makie.axislegend(ax, position=:lb, framevisible=false)

	ax = Makie.Axis(fig[2,1][1,3], xlabel="PC1", ylabel="PC2", width=250, height=250, xgridvisible=false, ygridvisible=false, title="Probed sequences")
	Makie.scatter!(ax, __proj_hits[:, end], __proj_hits[:, end - 1], markersize=10, color=(:gray, 0.5), label="MSA")
	Makie.scatter!(ax,
	    __proj_probed[(_probed_origin .== "RF00162_seed70") .| (_probed_origin .== "RF00162_full30"), end],
	    __proj_probed[(_probed_origin .== "RF00162_seed70") .| (_probed_origin .== "RF00162_full30"), end - 1],
	    markersize=10, color=:black, label="Natural", marker=:cross
	)
	Makie.scatter!(ax,
	    __proj_probed[_probed_origin .== "RF00162_syn_inf", end],
	    __proj_probed[_probed_origin .== "RF00162_syn_inf", end - 1],
	    markersize=10, color=:red, label="rCM", marker=:cross
	)
	Makie.scatter!(ax,
	    __proj_probed[_probed_origin .== "RF00162_syn_rbm", end],
	    __proj_probed[_probed_origin .== "RF00162_syn_rbm", end - 1],
	    markersize=10, color=:blue, label="RBM", marker=:cross
	)
	Makie.axislegend(ax, position=:lb, framevisible=false, patchlabelgap=-3)

	Makie.Label(fig[1,1][1,1][1,1,Makie.TopLeft()], "A)", font=:bold)
	Makie.Label(fig[1,1][1,2][1,1,Makie.TopLeft()], "B)", font=:bold)
	Makie.Label(fig[2,1][1,1][1,1,Makie.TopLeft()], "C)", font=:bold)
	Makie.Label(fig[2,1][1,2][1,1,Makie.TopLeft()], "D)", font=:bold)
	Makie.Label(fig[2,1][1,3][1,1,Makie.TopLeft()], "E)", font=:bold)

	# fig[0,4] = Makie.Label(fig, "Natural sequences", font=:bold)
	# fig[0,5] = Makie.Label(fig, "Generated sequences", font=:bold)

	Makie.resize_to_layout!(fig)
	#Makie.save("/workspaces/SamApp.jl/notebooks/2024-03-14 New paper figures/Figures/PCA.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═5d5eb95f-6801-47c3-a82e-afab35559c0e
# ╠═952a79aa-12fa-11ef-2326-d5fe00ff3dfe
# ╠═f01ccb33-26f8-4fb9-94a2-88be47db3762
# ╠═897ffcd7-e59a-4121-975d-29bc9a27988b
# ╠═49b0a181-8c3b-41dc-b3e0-9e55efa87bd9
# ╠═4a2a79bc-6fb1-45d3-ad4a-abfc889f2737
# ╠═5ca1d179-8dda-4e4b-9bab-b1cf4f6bafbe
# ╠═98c97d8b-9849-48dc-ae6c-c0924d14e560
# ╠═adbed31f-5777-4b6f-acf9-f37a472ef5d8
# ╠═de469256-6471-4fa9-854b-7b7783295aef
# ╠═a9fb22f9-4c47-4430-a0aa-cc4075b86f1f
# ╠═ebabc892-15f1-4b91-be5d-b85676d94fa7
# ╠═2a184e1f-bcba-48a3-a54b-793a25d9501e
# ╠═05b55f01-751c-4910-964f-aadb33f844c0
# ╠═a72c9391-c8f5-4bf5-94bd-3e9bbc326fa2
# ╠═bf22e1b2-f8bd-4990-938d-4db503660c3d
# ╠═ac2f96ec-b31f-4f63-8fd3-4390d687755b
# ╠═aa4f2963-ec25-4469-adea-f7f16c69db26
# ╠═aa8a002e-693a-402b-aaf1-576bf59e3a4e
# ╠═3a00d0fd-17e1-4518-8fca-ac1ccd0f7d28
# ╠═660fa41b-c6c2-4b41-a698-82530c8dc4ef
# ╠═8399bdac-bd99-4a20-be2c-2055db32a5a6
# ╠═fd9539ff-e725-4b18-94fd-e6a75fc2ec73
# ╠═6538d3e8-54d8-484d-8baa-3d3cc44ed72c
# ╠═71b73684-6053-4efc-92b7-391267f26ad9
# ╠═e4309e31-4086-4b3f-baec-8798ab975c44
# ╠═3a5e2efa-a739-42c0-a76c-50c6c07487e9
# ╠═0d4f2776-4666-48dd-80df-990c479569df
# ╠═488c249d-9810-4697-a6ca-7ab326b53762
# ╠═0eeacfc3-2240-46e2-b851-9d574771d1d7
# ╠═896d95be-2a0c-4c4b-9c73-d6e2d1452306
# ╠═0bd3ffae-c009-4740-aa57-4dd5126eb848
# ╠═7aba8c42-93a8-4528-90cd-b9a308960e80
# ╠═16508ca2-1c56-4136-842d-397ab4760586
# ╠═54d09468-5516-4568-9582-8172cc40c4f5
# ╠═1b1989e8-36d2-474a-aa6e-be02bebb8e00
# ╠═219bd907-7417-4718-aaf8-4555fb8f2e53
# ╠═cf5a9319-b99d-4c0a-a75c-f633fad52c63
# ╠═9f9e4969-bfbb-4f2e-971b-1265bf461a16
# ╠═92fec252-d3a6-429e-ba3d-01139abc7bc7
# ╠═66d621fa-231f-415b-aa0b-ea6dc88cbaf1
# ╠═68768649-9e5b-48e7-a33a-c50e20389a94
# ╠═51789884-f149-4562-9656-06b1e210fb43
# ╠═b0df20e8-cd4d-4800-a348-d7bbd3e3cd65
# ╠═93a22e9b-5f2f-4d1f-b693-e4489c6befa1
# ╠═703fde1e-33b2-4f66-b0b7-9e883f6470d4
# ╠═4e1d8cc4-e169-484a-9c08-ed56db0350f5
# ╠═ac689699-10a5-4989-bd27-46a0fae25be2
# ╠═373b6eb8-7e76-430a-8005-6eb688ba660e
# ╠═d056e7a9-3cc6-4179-9ce8-f8c567b1a56f
# ╠═fd502d5b-c54d-4f02-88bb-bfacdd6e350c
# ╠═d89ee7c9-f646-4832-bf32-427a9aa8a640
# ╠═38c55494-e38a-4358-bc30-ba0f2880eb10
# ╠═ea3ad48b-b0ac-4054-9548-fa9ed1842a16

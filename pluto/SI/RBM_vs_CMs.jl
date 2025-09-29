### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ f4e57f90-2c1f-49af-a246-0c6a765e1a12
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ b8692df8-26f5-4de4-9cbe-0a248cd55f59
using BioSequences: LongRNA

# ╔═╡ 50f1ea9c-c14e-4d8f-9048-c6647d8abee9
using Distributions: Gamma

# ╔═╡ ce363346-8eca-4e7d-9d61-47446db3cd67
using Makie: @L_str

# ╔═╡ 77892981-16dc-4402-b15a-3ab205ae93a4
using NaNStatistics: nanmean

# ╔═╡ e06158a1-f844-4e43-bfab-b8d4e0beac77
using NaNStatistics: nanstd

# ╔═╡ 52c1d889-080f-42b2-8cbc-a94be40ea43e
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2323766a-00fb-4226-bb57-77884fa3c43a
using Statistics: cor

# ╔═╡ dec068c4-eab2-4d05-bc9b-c19e7e8ee626
using Statistics: mean

# ╔═╡ 67082dd1-db8e-40e6-ae83-8029eb1221b1
using NaNStatistics: nansum

# ╔═╡ a9e38c8d-19ed-44a4-8374-a6bf8e857fb2
md"""
# Imports
"""

# ╔═╡ 7dd53826-d7c0-4bb3-93d4-20a3a87da167
import PlutoUI

# ╔═╡ dd79b460-6448-45c8-acb7-8cf0b7b03d3b
import CairoMakie, Makie

# ╔═╡ 6a1c2140-2b8b-4f58-ae4e-8f8dc473c089


# ╔═╡ 7bc08d72-e56c-4467-94a0-582dbded9ad4
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 86c7fcda-2d67-4b7f-a735-ce300c0227ff
import SamApp2025

# ╔═╡ b374aa45-88b5-413e-afb8-37cc97c9faaa
import StatsBase

# ╔═╡ 67d956b2-18c9-4e12-9d9e-e4a54b056f16
import Infernal

# ╔═╡ 72bb655b-e9db-4834-8d16-7aa51c6bba95
import Rfam

# ╔═╡ 7829c3b2-83d3-47bc-bbf7-e31751031bcd
import FASTX

# ╔═╡ bb893d2b-25aa-49ef-bd50-a2c6ada2448c
PlutoUI.TableOfContents()

# ╔═╡ 85c19ce6-4096-4331-aec2-082cc4c5f24a
md"""
# Load data
"""

# ╔═╡ 3e66fef9-612a-4a9b-89aa-eefb2b82cd8c
# RBM samples
sampled_v = SamApp2025.rbm2022samples(); # faster

# ╔═╡ 602fa5d2-88dd-48f7-ad67-f7331d63c473
# CM model from Rfam (this has the noisy floor!)
Rfam_cm = SamApp2025.rfam_RF00162_rfam_cm()

# ╔═╡ 952066b8-9430-43f4-8f27-895888875891
# these are already aligned and without inserts
RF00162_hits_sequences = SamApp2025.rfam_RF00162_hits();

# ╔═╡ 7bd2c5f8-16bf-4d98-aa19-96b97fc7332f
# aligned hits, used to train a new noiseless CM model (in Stockholm format, without inserts!)
RF00162_hits_stk = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true);

# ╔═╡ fbce70a0-3b12-4ce6-87e3-34bfbd884382
# fit new CM model using full alignment (without inserts), and without entropic noise
Denoised_cm = SamApp2025.rfam_RF00162_denoised_cm()

# ╔═╡ a89f2164-0baa-4613-bae9-941d04ef7a31
Rfam_cm_emitted_sequences = SamApp2025.infernal_cm_emit_sequences(Rfam_cm.out; N=5000, inserts=false)

# ╔═╡ 5d30684a-fb71-4621-9c2d-e7896dbfd50f
Denoised_cm_emitted_sequences = SamApp2025.infernal_cm_emit_sequences(Denoised_cm.cmout; N=5000, inserts=false)

# ╔═╡ 9ccdf111-1d2b-41e1-b330-5eb096f5f901
md"""
# Unknotted CM
"""

# ╔═╡ cd2bc167-d3dd-4bf8-ae33-608b22dcd9be
# fit new CM model using full alignment, with untangled pseudoknot (without inserts), and without entropic noise
Uknotted_cm_permutted = SamApp2025.rfam_RF00162_unknotted_cm()

# ╔═╡ 457541e2-2638-4fe1-a9ba-1b58d96da6c0
# Consensus secondary structure of RF00162 in WUSS format
wuss = SamApp2025.rfam_ss("RF00162")

# ╔═╡ f5d1965d-21b8-4395-ab37-6d59a7d5658f
perm = [
    1:findlast(==('A'), wuss); # first segment up to first branch of pseudoknot
    findfirst(==('a'), wuss):findlast(==('a'), wuss); # second branch of pseudoknot
    (findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1); # what's in between the pseudoknot
    (findlast(==('a'), wuss) + 1):108 # what's after the pseudoknot
];

# ╔═╡ 8dc56958-fac1-4971-ad6b-8f8c8c79fc48
# build training alignment for the Untangled CM
RF00162_hits_stk_permuted = tempname()

# ╔═╡ 01ba4b08-5135-48dc-bfa6-9cd7e8f9ee8e
# RF00162_hits_stk has the hits in Stockholm format, with one line per sequence
open(RF00162_hits_stk_permuted, "w") do file
    for (line_index, line) in enumerate(eachline(RF00162_hits_stk.out))
        if startswith(line, "# STOCKHOLM") || startswith(line, "#=GF") || startswith(line, "#=GS") || isempty(line) || line == "//"
            write(file, line, '\n') # these comment lines we just copy
        elseif startswith(line, "#=GR")
            continue # skip GR annotations
        elseif startswith(line, "#=GC SS_cons") # consensus secondary structure
            _ss = line[41:end]
            @assert _ss  == "((((((((,,,,<<<<<---<<<_____>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<----<<<<<<_____>>>>>>-->))))))))"
            @assert wuss == "((((((((,,,,<<<<<---<<<_AAAA>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<aaaa<<<<<<_____>>>>>>-->))))))))"
            @assert _ss == replace(wuss, 'A' => '_', 'a' => '-')
            # add parenthesis for pseudoknot
            _ss = _ss[1:(findfirst(==('A'), wuss) - 1)] * "((((" * _ss[(findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1)] * "))))" * _ss[(findlast(==('a'), wuss) + 1):end]
            @assert length(_ss) == 108
            write(file, line[1:40] * _ss[perm], '\n') # write untangled (permuted) secondary structure
        elseif startswith(line, "#=GC RF") # consensus sequence
            @assert line[41:end] == "cucUuAUcaAGAGgGGcgGAGGGAcuGGCCCuaUGAAgCCcCgGCAACCccccauaauaaggggAaGGUGCcAAuuCCugCcggccauuaaggccgGaaagAUaAgag"
            write(file, line[1:40] * line[41:end][perm], '\n')
        else # this is a sequence line
            @assert length(line[41:end]) == 108
            write(file, line[1:40] * line[41:end][perm], '\n')
        end
    end
end

# ╔═╡ 8a5fe11d-6857-47ca-921b-c6dd30d71dc3
Untangled_cm_permuted_emitted_sequences = SamApp2025.infernal_cm_emit_sequences(Uknotted_cm_permutted.cmout; N=5000, inserts=false)

# ╔═╡ 3f5b9040-6cb8-4445-a010-9215d8128456
# Permute back to correct column locations
Untangled_cm_emitted_sequences = [seq[invperm(perm)] for seq in Untangled_cm_permuted_emitted_sequences];

# ╔═╡ 9e2d19b1-581a-4529-8b25-6d537263224b
md"""
# Infernal scores
"""

# ╔═╡ 699d33dc-474b-4d15-b78c-5ab778c30a34
# Infernal scores of hits, using Rfam CM model
RF00162_hits_Rfam_cm_scores = SamApp2025.infernal_score_sequences(Rfam_cm.out, [replace(string(seq), '-' => "") for seq = RF00162_hits_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 1c94e498-2a12-4151-b252-8682026d7c22
# Infernal scores of hits, using Denoised CM model
RF00162_hits_Denoised_cm_scores = SamApp2025.infernal_score_sequences(Denoised_cm.cmout, [replace(string(seq), '-' => "") for seq = RF00162_hits_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 2584436c-018b-4769-9449-888868227e2b
# Infernal scores of hits, using Untangled CM model
RF00162_hits_Untangled_cm_scores = SamApp2025.infernal_score_sequences(Uknotted_cm_permutted.cmout, [replace(string(seq)[perm], '-' => "") for seq = RF00162_hits_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 66e0d54d-3c6c-4dd0-ac22-1c92d1357874
Rfam_cm_emitted_sequences_infernal_scores = SamApp2025.infernal_score_sequences(Rfam_cm.out, [replace(string(seq), '-' => "") for seq = Rfam_cm_emitted_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 077a471a-e69f-4ef9-af6b-61d76ca8f7f2
Denoised_cm_emitted_sequences_infernal_scores = SamApp2025.infernal_score_sequences(Denoised_cm.cmout, [replace(string(seq), '-' => "") for seq = Denoised_cm_emitted_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 788d547d-fc20-404f-940d-2fc977087ac9
Untangled_cm_emitted_sequences_infernal_scores = SamApp2025.infernal_score_sequences(Uknotted_cm_permutted.cmout, [replace(string(seq), '-' => "") for seq = Untangled_cm_permuted_emitted_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ f2ac783c-22ac-460a-a680-5be2455fa35c
RBM_samples_Rfam_CM_infernal_scores = SamApp2025.infernal_score_sequences(Rfam_cm.out, [replace(string(seq), '-' => "") for seq = SamApp2025.rnaseq(sampled_v)]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ c902284d-814d-42a9-aa51-f715a50a8934
RBM_samples_Denoised_CM_infernal_scores = SamApp2025.infernal_score_sequences(Denoised_cm.cmout, [replace(string(seq), '-' => "") for seq = SamApp2025.rnaseq(sampled_v)]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 3e35e9d4-efab-4c57-bc84-7888c92852d1
RBM_samples_Untangled_CM_infernal_scores = SamApp2025.infernal_score_sequences(Uknotted_cm_permutted.cmout, [replace(string(seq)[perm], '-' => "") for seq = SamApp2025.rnaseq(sampled_v)]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ e4a649b4-19ed-4573-9c2a-c122fd61e81f
md"""
# Figure
"""

# ╔═╡ 64c22582-6695-4057-bb93-ea331a24c8af
let fig = Makie.Figure()
	_sz = 350

	ax = Makie.Axis(fig[1,1], xlabel="Rfam CM score", ylabel="RBM score", width=_sz, height=_sz, xticks=-20:20:110, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Rfam_cm_scores, -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(RF00162_hits_sequences)), label="Natural", color=(:gray, 0.5), markersize=10)
	Makie.scatter!(ax,
	    Rfam_cm_emitted_sequences_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(Rfam_cm_emitted_sequences))[1:2000],
	    label="Rfam CM", color=:red, markersize=5)
	Makie.scatter!(ax,
	    RBM_samples_Rfam_CM_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), sampled_v)[1:2000],
	    label="RBM", color=:blue, markersize=5)
	Makie.xlims!(ax, -10, 110)
	Makie.axislegend(ax, position=:rb)

	ax = Makie.Axis(fig[1,2], xlabel="Denoised CM score", ylabel="RBM score", width=_sz, height=_sz, xticks=20:20:140, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Denoised_cm_scores, -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(RF00162_hits_sequences)), label="Natural", color=(:gray, 0.5), markersize=10)
	Makie.scatter!(ax,
	    Denoised_cm_emitted_sequences_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(Denoised_cm_emitted_sequences))[1:2000],
	    label="Denoised CM", color=:red, markersize=5)
	Makie.scatter!(ax,
	    RBM_samples_Denoised_CM_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), sampled_v)[1:2000],
	    label="RBM", color=:blue, markersize=5)
	Makie.xlims!(ax, 20, 140)
	Makie.axislegend(ax, position=:rb)

	ax = Makie.Axis(fig[1,3], xlabel="Unknotted CM score", ylabel="RBM score", width=_sz, height=_sz, xticks=20:20:140, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Untangled_cm_scores, -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(RF00162_hits_sequences)), label="Natural", color=(:gray, 0.5), markersize=10)
	Makie.scatter!(ax,
	    Untangled_cm_emitted_sequences_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), SamApp2025.onehot(Untangled_cm_emitted_sequences))[1:2000],
	    label="Unknotted CM", color=:red, markersize=5)
	Makie.scatter!(ax,
	    RBM_samples_Untangled_CM_infernal_scores[1:2000],
	    -RBMs.free_energy(SamApp2025.rbm2022(), sampled_v)[1:2000],
	    label="RBM", color=:blue, markersize=5)
	Makie.xlims!(ax, 20, 140)
	Makie.axislegend(ax, position=:rb)

	ax = Makie.Axis(fig[2,1], xlabel="Rfam CM score", ylabel="Denoised CM score", width=_sz, height=_sz, xticks=-20:20:110, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Rfam_cm_scores, RF00162_hits_Denoised_cm_scores, color=:black, markersize=5)
	ax = Makie.Axis(fig[2,2], xlabel="Rfam CM score", ylabel="Unknotted CM score", width=_sz, height=_sz, xticks=20:20:140, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Rfam_cm_scores, RF00162_hits_Untangled_cm_scores, color=:black, markersize=5)
	ax = Makie.Axis(fig[2,3], xlabel="Denoised CM score", ylabel="Unknotted CM score", width=_sz, height=_sz, xticks=20:20:140, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, RF00162_hits_Denoised_cm_scores, RF00162_hits_Untangled_cm_scores, color=:black, markersize=5)

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,2][1, 1, Makie.TopLeft()], "B)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,3][1, 1, Makie.TopLeft()], "C)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "D)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,2][1, 1, Makie.TopLeft()], "E)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,3][1, 1, Makie.TopLeft()], "F)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA/cossio/SAM/2024/SamApp2025.jl/pluto/SI/Figures/2024-04-15 RBMs vs. CMs.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═a9e38c8d-19ed-44a4-8374-a6bf8e857fb2
# ╠═f4e57f90-2c1f-49af-a246-0c6a765e1a12
# ╠═7dd53826-d7c0-4bb3-93d4-20a3a87da167
# ╠═dd79b460-6448-45c8-acb7-8cf0b7b03d3b
# ╠═6a1c2140-2b8b-4f58-ae4e-8f8dc473c089
# ╠═7bc08d72-e56c-4467-94a0-582dbded9ad4
# ╠═86c7fcda-2d67-4b7f-a735-ce300c0227ff
# ╠═b374aa45-88b5-413e-afb8-37cc97c9faaa
# ╠═67d956b2-18c9-4e12-9d9e-e4a54b056f16
# ╠═72bb655b-e9db-4834-8d16-7aa51c6bba95
# ╠═7829c3b2-83d3-47bc-bbf7-e31751031bcd
# ╠═b8692df8-26f5-4de4-9cbe-0a248cd55f59
# ╠═50f1ea9c-c14e-4d8f-9048-c6647d8abee9
# ╠═ce363346-8eca-4e7d-9d61-47446db3cd67
# ╠═77892981-16dc-4402-b15a-3ab205ae93a4
# ╠═e06158a1-f844-4e43-bfab-b8d4e0beac77
# ╠═52c1d889-080f-42b2-8cbc-a94be40ea43e
# ╠═2323766a-00fb-4226-bb57-77884fa3c43a
# ╠═dec068c4-eab2-4d05-bc9b-c19e7e8ee626
# ╠═67082dd1-db8e-40e6-ae83-8029eb1221b1
# ╠═bb893d2b-25aa-49ef-bd50-a2c6ada2448c
# ╠═85c19ce6-4096-4331-aec2-082cc4c5f24a
# ╠═3e66fef9-612a-4a9b-89aa-eefb2b82cd8c
# ╠═602fa5d2-88dd-48f7-ad67-f7331d63c473
# ╠═952066b8-9430-43f4-8f27-895888875891
# ╠═7bd2c5f8-16bf-4d98-aa19-96b97fc7332f
# ╠═fbce70a0-3b12-4ce6-87e3-34bfbd884382
# ╠═a89f2164-0baa-4613-bae9-941d04ef7a31
# ╠═5d30684a-fb71-4621-9c2d-e7896dbfd50f
# ╠═9ccdf111-1d2b-41e1-b330-5eb096f5f901
# ╠═cd2bc167-d3dd-4bf8-ae33-608b22dcd9be
# ╠═457541e2-2638-4fe1-a9ba-1b58d96da6c0
# ╠═f5d1965d-21b8-4395-ab37-6d59a7d5658f
# ╠═8dc56958-fac1-4971-ad6b-8f8c8c79fc48
# ╠═01ba4b08-5135-48dc-bfa6-9cd7e8f9ee8e
# ╠═8a5fe11d-6857-47ca-921b-c6dd30d71dc3
# ╠═3f5b9040-6cb8-4445-a010-9215d8128456
# ╠═9e2d19b1-581a-4529-8b25-6d537263224b
# ╠═699d33dc-474b-4d15-b78c-5ab778c30a34
# ╠═1c94e498-2a12-4151-b252-8682026d7c22
# ╠═2584436c-018b-4769-9449-888868227e2b
# ╠═66e0d54d-3c6c-4dd0-ac22-1c92d1357874
# ╠═077a471a-e69f-4ef9-af6b-61d76ca8f7f2
# ╠═788d547d-fc20-404f-940d-2fc977087ac9
# ╠═f2ac783c-22ac-460a-a680-5be2455fa35c
# ╠═c902284d-814d-42a9-aa51-f715a50a8934
# ╠═3e35e9d4-efab-4c57-bc84-7888c92852d1
# ╠═e4a649b4-19ed-4573-9c2a-c122fd61e81f
# ╠═64c22582-6695-4057-bb93-ea331a24c8af

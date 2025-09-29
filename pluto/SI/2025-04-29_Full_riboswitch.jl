### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 9e5ab2f6-2505-11f0-2a08-2599aa3bd022
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 75226411-330d-4450-a6d2-0ca3d6774f65
using Statistics: mean, cov, std, var

# ╔═╡ 89c7b400-aa05-4a38-8dbb-130eec789084
using Distributions: PoissonBinomial, pdf, logpdf

# ╔═╡ b97d20b6-70e2-4194-b49c-3166761d6403
using LogExpFunctions: xlogx

# ╔═╡ 73398a3e-1c08-4416-ac15-9bbe282834ce
using LinearAlgebra: diagind

# ╔═╡ 48c28fd4-1b35-4632-9215-a24f3ecb1a80
using Unitful: ustrip

# ╔═╡ c89c23fe-3b38-415b-b11d-541e814d94fb
md"# Imports"

# ╔═╡ 81a424bf-f09b-480a-99b3-528a561ecfad
import Makie, CairoMakie, Logomaker

# ╔═╡ f3ec55d7-80e5-467d-8aa6-7f7804dd5f58
import ViennaRNA, BioSequences, SamApp2025

# ╔═╡ ccd05903-4194-4d11-9025-1aac33df4fde
md"# Data"

# ╔═╡ 4083a737-315c-467c-9b0a-eb6e829e47dd
fullseqs = SamApp2025.artifact_load_20250428_aligned_full_riboswitch_sequences()

# ╔═╡ 6bd76c54-45cc-4158-b6c0-604503a2fb28
alnseqs = map(BioSequences.LongRNA{4}, [filter(!islowercase, replace(seq, '.' => "")) for seq = fullseqs])

# ╔═╡ 7f6c9a0b-7bf4-4e07-8226-38de3746420a
X = SamApp2025.onehot(map(BioSequences.LongRNA{4}, alnseqs));

# ╔═╡ fafcc597-a322-48a2-ad6d-eb9ef95b6527
wuss = SamApp2025.rfam_ss("RF00162"; inserts=false)

# ╔═╡ 84cb87f0-e81a-4799-aa19-57b9b5fb63c3
ss = SamApp2025.clean_wuss(wuss)

# ╔═╡ f1de1a29-8301-447f-adad-ec906d94724a
ss_without_P1 = join([i ∈ SamApp2025.RF00162_sites_annotated_secondary_structure().p1 ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ a0ee4304-50ff-440a-b363-2ea91a911384
# full seq, with P1 and terminator bound (OFF state)
ss_full_P1_and_terminator = ss * repeat('.', length(109:115)) * repeat('(', length(116:123)) * repeat('.', length(124:130)) * repeat(')', length(131:138)) * repeat('.', length(139:153))

# ╔═╡ 85895b80-369c-474e-a134-9b95b45b60a6
# full seq, with P1 but without terminator
ss_full_P1_and_no_terminator = ss * repeat('.', length(109:153))

# ╔═╡ d1244e1c-91f8-4378-8734-c52d3fa95385
# full seq, with P1 unbound and anti-terminator bound (ON state)
ss_full_noP1 = repeat('.', 8) * ss[9:100] * repeat('(', 8) *
	repeat('.', length(109:115)) * repeat(')', length(116:123)) * repeat('.', length(124:130)) * repeat('.', length(131:138)) * repeat('.', length(139:153))

# ╔═╡ a6744f6a-477d-47c4-9dcf-6466026893a5
Vienna_energies_aptamer_only_ss = [ustrip(ViennaRNA.energy(string(seq[1:108]), ss)) for seq = alnseqs];

# ╔═╡ 93c341ff-dde3-4dd9-b245-fc9db2d9a5ab
Vienna_energies_aptamer_only_noP1 = [ustrip(ViennaRNA.energy(string(seq[1:108]), ss_without_P1)) for seq = alnseqs];

# ╔═╡ 5eefdb8e-85bd-4116-9346-e9394de89475
Vienna_energies_full_P1_and_terminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_P1_and_terminator)) for seq = alnseqs];

# ╔═╡ 06236adb-22da-4e8e-80f6-f4398e1f5b5f
Vienna_energies_full_P1_and_no_terminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_P1_and_no_terminator)) for seq = alnseqs];

# ╔═╡ e66a8b94-73a8-42a4-b65f-4190d37fe172
Vienna_energies_full_noP1_antiterminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_noP1)) for seq = alnseqs];

# ╔═╡ c5735f3f-8982-4eab-bc30-39b0b80a27c5
P1_deltaF_aptamer_only = Vienna_energies_aptamer_only_ss - Vienna_energies_aptamer_only_noP1

# ╔═╡ 4a8133b7-f404-4364-ae3a-ab7457364456
P1_deltaF_full = Vienna_energies_full_P1_and_terminator - Vienna_energies_full_noP1_antiterminator

# ╔═╡ eddb9d24-017e-40a6-ae17-00cc4e407e00
P1_deltaF_full_no_terminator = Vienna_energies_full_P1_and_no_terminator - Vienna_energies_full_noP1_antiterminator

# ╔═╡ 0fbc4fa8-5903-4afb-9759-2ac24a637003
BPs_full = [ViennaRNA.bpp(string(seq)) for seq = alnseqs];

# ╔═╡ dd32923e-e86f-4165-8d38-5a7cdfb75bd6
BPs_aptamer = [ViennaRNA.bpp(string(seq[1:108])) for seq = alnseqs];!

# ╔═╡ 43564b0d-3a29-4f88-8d15-8120300660fa
let seq = string(alnseqs[30])
	println(seq[1:108])
	println(ViennaRNA.mfe(seq[1:108])[1])
	println(seq)
	println(ViennaRNA.mfe(seq)[1])
end

# ╔═╡ 888345b7-8661-47d5-bedd-5d59612071b6
let fig = Makie.Figure()
	bins = 0:0.1:1
	for i = 1:8
		j = 109 - i
		ax = Makie.Axis(fig[fld1(i, 4), mod1(i, 4)]; width=150, height=150, xlabel="Base-pair prob. P1", ylabel="aptamers (freq.)", title="pair = ($i, $j)")
		Makie.hist!(ax, [bp[i,j] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i,j] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)
		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 10)
	end

	_xs = 0:8
	bp_probs_aptamer = [pdf(PoissonBinomial([bp[i, 109 - i] for i=1:8]), n) for n = _xs, bp = BPs_aptamer]
	bp_probs_riboswitch = [pdf(PoissonBinomial([bp[i, 109 - i] for i=1:8]), n) for n = _xs, bp = BPs_full]

	ax = Makie.Axis(fig[1,5], width=150, height=150, xlabel="Num. paired sites", ylabel="Average probability")
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)), vec(std(bp_probs_aptamer; dims=2))/2; color=:blue, linewidth=1, whiskerwidth=5)
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_riboswitch; dims=2)), vec(std(bp_probs_riboswitch; dims=2))/2; color=:black, linewidth=1, whiskerwidth=5)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)); label="aptamer only", color=:blue, linewidth=2)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_riboswitch; dims=2)); label="full riboswitch", color=:black, linewidth=2)
	Makie.ylims!(ax, 0, 0.5)

	ax = Makie.Axis(fig[2,5], width=150, height=150, xlabel="Most likely no. of paired sites", ylabel="aptamers (freq.)")
	Makie.stephist!(ax, [argmax(logpdf(PoissonBinomial([bp[i, 109 - i] for i=1:8]), n) for n = 0:10) - 1 for bp = BPs_aptamer]; normalization=:pdf, label="aptamer", color=:blue, linewidth=2, bins=-0.5:1:8.5)
	Makie.stephist!(ax, [argmax(logpdf(PoissonBinomial([bp[i, 109 - i] for i=1:8]), n) for n = 0:10) - 1 for bp = BPs_full]; normalization=:pdf, label="riboswitch", color=:black, linewidth=2, bins=-0.5:1:8.5)
	Makie.ylims!(ax, 0, 0.5)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	# ax = Makie.Axis(fig[1:2, 5], width=350, height=350, xlabel="Prob. all unpaired (log10)", ylabel="aptamers (freq.)")
	# Makie.hist!(ax, [log10(prod(1-bp[i, 109 - i] for i=1:8)) for bp = BPs_full]; normalization=:pdf, label="full", color=:gray, bins=-23:2:0)
	# Makie.stephist!(ax, [log10(prod(1-bp[i, 109 - i] for i=1:8)) for bp = BPs_aptamer]; normalization=:pdf, label="aptamer only", color=:blue, linewidth=2, bins=-23:2:0)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i, 109 - i]) for i=1:8) for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i, 109 - i]) for i=1:8) for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)
	# Makie.ylims!(ax, 0, 0.15)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 858cf27e-ea29-4d4f-9205-c6de4ba49f33
let fig = Makie.Figure()
	bins = 0:0.1:1

	is = [13:17; 21:23]
	js = reverse([29:31; 38:42])

	for (n, (i,j)) = enumerate(zip(is, js))
		ax = Makie.Axis(fig[fld1(n, 4), mod1(n, 4)]; width=150, height=150, xlabel="Base-pair prob. (P2)", ylabel="aptamers (freq.)", title="pair = ($i, $j)")

		Makie.hist!(ax, [bp[i,j] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i,j] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)

		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)

		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 7)
	end

	_xs = 0:8
	bp_probs_aptamer = [pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = _xs, bp = BPs_aptamer]
	bp_probs_full = [pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = _xs, bp = BPs_full]

	ax = Makie.Axis(fig[1,5], width=150, height=150, xlabel="Num. paired sites", ylabel="Average probability")
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)), vec(std(bp_probs_aptamer; dims=2))/2; color=:blue, linewidth=1, whiskerwidth=5)
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_full; dims=2)), vec(std(bp_probs_full; dims=2))/2; color=:black, linewidth=1, whiskerwidth=5)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)); label="aptamer only", color=:blue, linewidth=2)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_full; dims=2)); label="full riboswitch", color=:black, linewidth=2)
	Makie.ylims!(ax, 0, 0.5)

	ax = Makie.Axis(fig[2,5], width=150, height=150, xlabel="Most likely no. of paired sites", ylabel="aptamers (freq.)")
	Makie.stephist!(ax, [argmax(pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = 0:10) - 1 for bp = BPs_aptamer]; normalization=:pdf, label="aptamer", color=:blue, linewidth=2, bins=-0.5:1:8.5)
	Makie.stephist!(ax, [argmax(pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = 0:10) - 1 for bp = BPs_full]; normalization=:pdf, label="riboswitch", color=:black, linewidth=2, bins=-0.5:1:8.5)
	Makie.ylims!(ax, 0, 0.5)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	# ax = Makie.Axis(fig[1:2, 5], width=350, height=350, xlabel="Prob. of any unpaired site (log10)", ylabel="aptamers (freq.)")
	# Makie.hist!(ax, [log10(prod(1-bp[i,j] for (i,j)=zip(is,js))) for bp = BPs_full]; normalization=:pdf, label="full", color=:gray, bins=-23:2:0)
	# Makie.stephist!(ax, [log10(prod(1-bp[i,j] for (i,j)=zip(is,js))) for bp = BPs_aptamer]; normalization=:pdf, label="aptamer only", color=:blue, linewidth=2, bins=-23:2:0)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i,j]) for (i,j)=zip(is,js)) for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i,j]) for (i,j)=zip(is,js)) for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)
	# Makie.ylims!(ax, 0, 0.15)
	# Makie.axislegend(ax; position=:lt)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ df7fccd6-a5c4-47c4-afe4-a55a19ed57e3
let fig = Makie.Figure()
	bins = 0:0.1:1
	is = 81:86
	js = 97 .- is .+ 81
	for (i,j) = zip(is,js)
		n = i - 80
		ax = Makie.Axis(fig[fld1(n,3), mod1(n,3)]; width=150, height=150, xlabel="Base-pair prob. (P4)", ylabel="aptamers (freq.)", title="pair = ($i, $j)")

		Makie.hist!(ax, [bp[i,j] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i,j] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)

		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)

		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 10)
	end

	_xs = 0:6
	bp_probs_aptamer = [pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = _xs, bp = BPs_aptamer]
	bp_probs_full = [pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = _xs, bp = BPs_full]

	ax = Makie.Axis(fig[1,4], width=150, height=150, xlabel="Num. paired sites", ylabel="Average probability")
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)), vec(std(bp_probs_aptamer; dims=2))/2; color=:blue, linewidth=1, whiskerwidth=5)
	Makie.errorbars!(ax, _xs, vec(mean(bp_probs_full; dims=2)), vec(std(bp_probs_full; dims=2))/2; color=:black, linewidth=1, whiskerwidth=5)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_aptamer; dims=2)); label="aptamer only", color=:blue, linewidth=2)
	Makie.scatterlines!(ax, _xs, vec(mean(bp_probs_full; dims=2)); label="full riboswitch", color=:black, linewidth=2)
	Makie.ylims!(ax, 0, 0.5)

	ax = Makie.Axis(fig[2,4], width=150, height=150, xlabel="Most likely no. of paired sites", ylabel="aptamers (freq.)")
	Makie.stephist!(ax, [argmax(pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = 0:7) - 1 for bp = BPs_aptamer]; normalization=:pdf, label="aptamer", color=:blue, linewidth=2, bins=-0.5:1:6.5)
	Makie.stephist!(ax, [argmax(pdf(PoissonBinomial([bp[i,j] for (i,j)=zip(is,js)]), n) for n = 0:7) - 1 for bp = BPs_full]; normalization=:pdf, label="riboswitch", color=:black, linewidth=2, bins=-0.5:1:6.5)
	Makie.ylims!(ax, 0, 0.5)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	# ax = Makie.Axis(fig[1:2, 4], width=350, height=350, xlabel="Prob. of any unpaired site (log10)", ylabel="aptamers (freq.)")
	# Makie.hist!(ax, [log10(prod(1-bp[i,97-i+81] for i=81:86)) for bp = BPs_full]; normalization=:pdf, label="full", color=:gray, bins=-23:2:0)
	# Makie.stephist!(ax, [log10(prod(1-bp[i,97-i+81] for i=81:86)) for bp = BPs_aptamer]; normalization=:pdf, label="aptamer only", color=:blue, linewidth=2, bins=-23:2:0)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i,97-i+81]) for i=81:86) for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
	# Makie.vlines!(ax, mean([sum(log10(1-bp[i,97-i+81]) for i=81:86) for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)
	# Makie.ylims!(ax, 0, 0.15)
	# Makie.axislegend(ax; position=:lt)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═c89c23fe-3b38-415b-b11d-541e814d94fb
# ╠═9e5ab2f6-2505-11f0-2a08-2599aa3bd022
# ╠═81a424bf-f09b-480a-99b3-528a561ecfad
# ╠═f3ec55d7-80e5-467d-8aa6-7f7804dd5f58
# ╠═75226411-330d-4450-a6d2-0ca3d6774f65
# ╠═89c7b400-aa05-4a38-8dbb-130eec789084
# ╠═b97d20b6-70e2-4194-b49c-3166761d6403
# ╠═73398a3e-1c08-4416-ac15-9bbe282834ce
# ╠═48c28fd4-1b35-4632-9215-a24f3ecb1a80
# ╠═ccd05903-4194-4d11-9025-1aac33df4fde
# ╠═4083a737-315c-467c-9b0a-eb6e829e47dd
# ╠═6bd76c54-45cc-4158-b6c0-604503a2fb28
# ╠═7f6c9a0b-7bf4-4e07-8226-38de3746420a
# ╠═fafcc597-a322-48a2-ad6d-eb9ef95b6527
# ╠═84cb87f0-e81a-4799-aa19-57b9b5fb63c3
# ╠═f1de1a29-8301-447f-adad-ec906d94724a
# ╠═a0ee4304-50ff-440a-b363-2ea91a911384
# ╠═85895b80-369c-474e-a134-9b95b45b60a6
# ╠═d1244e1c-91f8-4378-8734-c52d3fa95385
# ╠═a6744f6a-477d-47c4-9dcf-6466026893a5
# ╠═93c341ff-dde3-4dd9-b245-fc9db2d9a5ab
# ╠═5eefdb8e-85bd-4116-9346-e9394de89475
# ╠═06236adb-22da-4e8e-80f6-f4398e1f5b5f
# ╠═e66a8b94-73a8-42a4-b65f-4190d37fe172
# ╠═c5735f3f-8982-4eab-bc30-39b0b80a27c5
# ╠═4a8133b7-f404-4364-ae3a-ab7457364456
# ╠═eddb9d24-017e-40a6-ae17-00cc4e407e00
# ╠═0fbc4fa8-5903-4afb-9759-2ac24a637003
# ╠═dd32923e-e86f-4165-8d38-5a7cdfb75bd6
# ╠═43564b0d-3a29-4f88-8d15-8120300660fa
# ╠═888345b7-8661-47d5-bedd-5d59612071b6
# ╠═858cf27e-ea29-4d4f-9205-c6de4ba49f33
# ╠═df7fccd6-a5c4-47c4-afe4-a55a19ed57e3

### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ cb9ef902-4a25-41f1-bd09-f9bb7f5297e0
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ c65a6736-34b8-46f9-b72f-0620a7d9b802
using BioSequences: LongRNA

# ╔═╡ 6dcf3a81-2081-4a7f-9294-ef8199c0ea98
using Statistics: mean

# ╔═╡ 2e8b3992-a298-4fdf-a499-e33608e844f4
using DataFrames: DataFrame

# ╔═╡ 3dda1c7a-3bb2-45d5-b37e-f2cfa2aaa24b
using LogExpFunctions: xlogx

# ╔═╡ eee444c1-f6ba-4aa4-bb98-22b5a8d068ca
md"# Imports"

# ╔═╡ 7227bb0f-77a3-4707-8c49-e742b3e0a1d8
import Makie

# ╔═╡ 6506a7eb-e4cb-471b-b0a5-ebd76a8b7262
import CairoMakie

# ╔═╡ 1e6598aa-0b5f-4143-9c95-991deaeaaab1
import FASTX

# ╔═╡ b6a4a192-9751-48d9-8b10-42b2e23755bc
import Infernal

# ╔═╡ dcb5958e-7a9d-4df7-82c5-145019736c7d
import SamApp2025

# ╔═╡ b7ea3497-5307-4728-88d6-f7afc02c26a8
import CSV

# ╔═╡ c91b17d8-7f6e-463e-b79b-6f554697878e
import Rfam

# ╔═╡ fa739257-4193-4e6c-b950-db161d5f62bd
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 49b52aea-cc64-4a40-a26f-b7e6c5c03527
import PlutoUI

# ╔═╡ 88b56aee-8fd0-471a-a643-c99005125ab8
import Logomaker, PythonPlot

# ╔═╡ 0d0fce23-94dc-4a5d-9bab-38f25c78e516
PlutoUI.TableOfContents()

# ╔═╡ 8d920fe4-3402-47a3-9943-e8b4ed3b256d
md"# Load data"

# ╔═╡ 1209ca4c-7f0a-4ea0-8170-6295b8b18ab4
msa = SamApp2025.rfam_RF00162_hits()

# ╔═╡ d581c996-4d93-4fed-8838-3487010bb715
md"# Sequence logo"

# ╔═╡ 9f701d85-3410-4049-8582-665451746c10
xlog2x(x) = xlogx(x) / log(oftype(x,2))

# ╔═╡ 275cd9fe-9f2c-486a-b6f2-b5077ff9edd0
NTs = collect("ACGU⊟")

# ╔═╡ ce46177d-d4df-433b-87ab-b96a76b66b25
seqlogo_color_scheme = Logomaker.color_scheme(
	'C' => "blue",
	'U' => "red",
	'A' => "green",
	'G' => "orange",
	'⊟' => "black"
)


# ╔═╡ b6a9b93e-33fe-4bbb-b0bf-025c72b33bd3
function seqlogo_entropic(p::AbstractMatrix; max_ylim=true)
    @assert size(p, 1) == 5 # nucleotides + gap
    w = p ./ sum(p; dims=1)
    H = sum(-xlog2x.(w); dims=1)
    @assert all(0 .≤ H .≤ log2(5))

    cons = w .* (log2(5) .- H)

	fig, ax = PythonPlot.subplots(1, 1, figsize=[15, 4])
	ax = PythonPlot.subplot(1, 1, 1)

    logo = Logomaker.Logo(cons, collect("ACGU⊟"); color_scheme=seqlogo_color_scheme, ax=ax)
    max_ylim && logo.ax.set_ylim(0, log2(5))
    logo.ax.set_ylabel("conservation (bits)")
    logo.ax.set_xlabel("site")
	logo.ax.set_xticks(5:5:108)
	logo.fig.tight_layout()
	logo.fig
    return logo
end

# ╔═╡ e5e3700b-f9a7-46cc-b5ef-e861837e8143
seqlogo_entropic(reshape(mean(SamApp2025.rbm2022samples(); dims=3), 5, 108)).fig

# ╔═╡ 82d312e7-4a3c-4594-a26b-317a7c64b1c2
seqlogo_entropic(reshape(mean(SamApp2025.onehot(msa); dims=3), 5, 108)).fig

# ╔═╡ Cell order:
# ╠═eee444c1-f6ba-4aa4-bb98-22b5a8d068ca
# ╠═cb9ef902-4a25-41f1-bd09-f9bb7f5297e0
# ╠═7227bb0f-77a3-4707-8c49-e742b3e0a1d8
# ╠═6506a7eb-e4cb-471b-b0a5-ebd76a8b7262
# ╠═1e6598aa-0b5f-4143-9c95-991deaeaaab1
# ╠═b6a4a192-9751-48d9-8b10-42b2e23755bc
# ╠═dcb5958e-7a9d-4df7-82c5-145019736c7d
# ╠═b7ea3497-5307-4728-88d6-f7afc02c26a8
# ╠═c91b17d8-7f6e-463e-b79b-6f554697878e
# ╠═fa739257-4193-4e6c-b950-db161d5f62bd
# ╠═49b52aea-cc64-4a40-a26f-b7e6c5c03527
# ╠═88b56aee-8fd0-471a-a643-c99005125ab8
# ╠═c65a6736-34b8-46f9-b72f-0620a7d9b802
# ╠═6dcf3a81-2081-4a7f-9294-ef8199c0ea98
# ╠═2e8b3992-a298-4fdf-a499-e33608e844f4
# ╠═3dda1c7a-3bb2-45d5-b37e-f2cfa2aaa24b
# ╠═0d0fce23-94dc-4a5d-9bab-38f25c78e516
# ╠═8d920fe4-3402-47a3-9943-e8b4ed3b256d
# ╠═1209ca4c-7f0a-4ea0-8170-6295b8b18ab4
# ╠═d581c996-4d93-4fed-8838-3487010bb715
# ╠═9f701d85-3410-4049-8582-665451746c10
# ╠═275cd9fe-9f2c-486a-b6f2-b5077ff9edd0
# ╠═ce46177d-d4df-433b-87ab-b96a76b66b25
# ╠═b6a9b93e-33fe-4bbb-b0bf-025c72b33bd3
# ╠═e5e3700b-f9a7-46cc-b5ef-e861837e8143
# ╠═82d312e7-4a3c-4594-a26b-317a7c64b1c2

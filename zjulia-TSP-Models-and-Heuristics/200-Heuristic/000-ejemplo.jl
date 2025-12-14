### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 36184f24-de30-11ec-3a6f-6962a488690c
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    Pkg.add([
        Pkg.PackageSpec("Plots")
        Pkg.PackageSpec("PlutoUI")
    ])
    Pkg.status()
end

# ╔═╡ 9260ef64-e782-431e-80a2-fe64ae002fe6
using Plots

# ╔═╡ ab17dac6-4197-4320-ac38-16072550731b
using Random

# ╔═╡ 4b030fb4-5bd0-47eb-9c9d-2815b595acab
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ 3cc260c6-dc49-4939-9fd1-0e616556d167
include(joinpath(@__DIR__, "100-function.jl"))

# ╔═╡ 2348b289-9f91-4e67-a365-5f346f01db01
include(joinpath(@__DIR__, "200-initial.jl"))

# ╔═╡ 1acdf778-7b93-493b-b6a6-5230bf9cda52
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 4d8f5b53-7607-48f5-95f5-ce6ba2dfaa0c
n = 7

# ╔═╡ cc282f51-0086-413f-94d9-2d4a5bba34e0
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ 222a9beb-0bf8-459d-9df1-115196ec90b3
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ 31b71c37-8bd1-4ce4-8789-e13f7705122d
# tour = tspRND(dist)

# tour = tspNN(dist)
# tour = tspNN(dist, s = 1)

# tour = tspCI(dist)
# tour = tspCI(dist, s = 1)

# tour = tspCWs(dist)
# tour = tspCWs(dist, lamda = 1.0)
# tour = tspCWs(dist, s = 1)
# tour = tspCWs(dist, lamda = 1.0, s = 1)

# tour = tspCWp(dist)
# tour = tspCWp(dist, lamda = 1.0)
# tour = tspCWp(dist, s = 1)
# tour = tspCWp(dist, lamda = 1.0, s = 1)

tour = tspNN(dist, s = 1)

# ╔═╡ 1b9e0fce-6d96-4f2f-8d98-ae5af461e98c
z = tspDist(dist, tour)

# ╔═╡ c3968752-aa12-4a7f-8a04-844259a4d85b
p = tspPlot(X, Y, tour)

# ╔═╡ Cell order:
# ╠═36184f24-de30-11ec-3a6f-6962a488690c
# ╠═1acdf778-7b93-493b-b6a6-5230bf9cda52
# ╠═9260ef64-e782-431e-80a2-fe64ae002fe6
# ╠═ab17dac6-4197-4320-ac38-16072550731b
# ╠═3cc260c6-dc49-4939-9fd1-0e616556d167
# ╠═2348b289-9f91-4e67-a365-5f346f01db01
# ╠═4d8f5b53-7607-48f5-95f5-ce6ba2dfaa0c
# ╠═cc282f51-0086-413f-94d9-2d4a5bba34e0
# ╠═222a9beb-0bf8-459d-9df1-115196ec90b3
# ╠═31b71c37-8bd1-4ce4-8789-e13f7705122d
# ╠═1b9e0fce-6d96-4f2f-8d98-ae5af461e98c
# ╠═c3968752-aa12-4a7f-8a04-844259a4d85b
# ╠═4b030fb4-5bd0-47eb-9c9d-2815b595acab

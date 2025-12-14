### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 51bfede6-1436-11ec-0cc6-bfff507e0e99
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

# ╔═╡ 879c42ee-f980-42e9-bb81-2f27672aecbf
using Plots

# ╔═╡ 5724ac80-d48f-434d-9c57-51bb2d655ea9
using Random

# ╔═╡ 91f30def-bd7d-4064-92d2-52f42c1159c4
using PlutoUI

# ╔═╡ 11fb6a38-ce89-41c5-bf72-477eec5d7b37
include(joinpath(@__DIR__, "../200-Heuristic/100-function.jl"))

# ╔═╡ fa0616ad-8697-4948-8db0-e5bc52309e84
include(joinpath(@__DIR__, "../200-Heuristic/200-initial.jl"))

# ╔═╡ 6c1d2ef0-cf3d-4993-88b8-6c000ba20ae3
include(joinpath(@__DIR__, "../300-LocalSearch/100-neighborhood.jl"))

# ╔═╡ da1b41b8-e924-423d-b76e-32a715b1f178
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 151335ef-2181-484b-a47d-b5948d3acb5e
n = 7

# ╔═╡ 26f6beff-c080-4294-b4b5-65a243fcefa7
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ 4c60dcfc-cfe8-4c09-a0f3-62674e59f8bb
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ 082943ba-3442-475b-b585-e913bfa34bb0
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

tour = tspRND(dist)

# ╔═╡ d049dd39-4736-4ff0-b2b0-981c61ee3296
md"""
## Simulated Annealing
"""

# ╔═╡ fab77d3b-203d-4033-bdd0-6ee1ff3cc572
function mhSA(tour::Array{Int,1}, dist::Array{Float64,2}, nbh::Array{Function,1})
    n, _ = size(dist)
    ctour = copy(tour)
    besttour = copy(tour)

    temp = 1_000.0
    temp_alpha = 0.9753
    maxiter = 50 * n

    busqueda0 = [tspDist(dist, ctour)] # actual ztour
    busqueda1 = copy(busqueda0)        # mejor  ztour
    busquedat = [temp]                 # actual temp
    busquedafn = ["base"]              # actual fn

    for _ in 1:maxiter
        nh = rand(nbh)
        newtour = nh(ctour, dist)
        newtourz = tspDist(dist, newtour)

        delta = newtourz - busqueda0[end]
        acepta = delta < 0.0 ? true : rand() < exp(-delta / temp)
        if acepta
            ctour = newtour
            push!(busqueda0, newtourz)
        else
            push!(busqueda0, busqueda0[end])
        end

        if newtourz < busqueda1[end]
            besttour = newtour
            push!(busqueda1, newtourz)
        else
            push!(busqueda1, busqueda1[end])
        end

        temp = temp_alpha * temp
        push!(busquedat, temp)
        push!(busquedafn, string(nh))
        @show busqueda0[end], nh
    end

    plot(busqueda0)
    plot!(busqueda1)

    savefig(joinpath(@__DIR__, "tsp-plot-SAz.png"))
    savefig(plot(busquedat), joinpath(@__DIR__, "tsp-plot-SAt.png"))

    return besttour
end

# ╔═╡ 526a7132-103b-4c0b-9bf5-7125de50320e
# nbh::Array{Function,1} = [tspInsSA]
# nbh::Array{Function,1} = [tspSwapSA]
# nbh::Array{Function,1} = [tsp2OptSA]
# nbh::Array{Function,1} = [tspInsSA, tspSwapSA, tsp2OptSA]
# nbh::Array{Function,1} = [tsp2OptSA, tspSwapSA, tspInsSA]
nbh::Array{Function,1} = [tspSwapSA, tspInsSA, tsp2OptSA]

# ╔═╡ fe14fe38-c519-4103-bf59-207cfd2eb402
tourMH = mhSA(tour, dist, nbh)

# ╔═╡ f5050494-1975-4558-9a1f-2a44b388b3dd
tspPlot(X, Y, tourMH)

# ╔═╡ 1e099701-6fd7-42d3-a76e-0077c4c8a174
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-SAz.png"))

# ╔═╡ 0c85af3d-16bc-493c-a4e2-43af9a61026c
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-SAt.png"))

# ╔═╡ 357e9cc4-a7bf-4618-8f1f-fc7af0601d7a
begin
    # using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ Cell order:
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═da1b41b8-e924-423d-b76e-32a715b1f178
# ╠═879c42ee-f980-42e9-bb81-2f27672aecbf
# ╠═5724ac80-d48f-434d-9c57-51bb2d655ea9
# ╠═91f30def-bd7d-4064-92d2-52f42c1159c4
# ╠═11fb6a38-ce89-41c5-bf72-477eec5d7b37
# ╠═fa0616ad-8697-4948-8db0-e5bc52309e84
# ╠═6c1d2ef0-cf3d-4993-88b8-6c000ba20ae3
# ╠═151335ef-2181-484b-a47d-b5948d3acb5e
# ╠═26f6beff-c080-4294-b4b5-65a243fcefa7
# ╠═4c60dcfc-cfe8-4c09-a0f3-62674e59f8bb
# ╠═082943ba-3442-475b-b585-e913bfa34bb0
# ╠═d049dd39-4736-4ff0-b2b0-981c61ee3296
# ╠═fab77d3b-203d-4033-bdd0-6ee1ff3cc572
# ╠═526a7132-103b-4c0b-9bf5-7125de50320e
# ╠═fe14fe38-c519-4103-bf59-207cfd2eb402
# ╠═f5050494-1975-4558-9a1f-2a44b388b3dd
# ╠═1e099701-6fd7-42d3-a76e-0077c4c8a174
# ╠═0c85af3d-16bc-493c-a4e2-43af9a61026c
# ╠═357e9cc4-a7bf-4618-8f1f-fc7af0601d7a

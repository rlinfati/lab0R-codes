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

# ╔═╡ 03fcd07b-99d5-4e2d-bc7f-9e00909139c6
using Random

# ╔═╡ 23121b1a-e4b1-4b47-8a2a-05fb11640818
using PlutoUI

# ╔═╡ e03fba74-b7f1-4203-a8a2-faa795a514b8
include(joinpath(@__DIR__, "../200-Heuristic/100-function.jl"))

# ╔═╡ da135245-37a5-40a1-b29b-2ba403f3af47
include(joinpath(@__DIR__, "../200-Heuristic/200-initial.jl"))

# ╔═╡ 28e27c70-5d42-4096-a40d-ba64848865e4
include(joinpath(@__DIR__, "../300-LocalSearch/100-neighborhood.jl"))

# ╔═╡ 5329815c-850b-4d8b-9846-a4ffdf02b448
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 3817cb46-09bb-41ff-8dd1-d345bc8f5136
n = 7

# ╔═╡ a30695c5-2057-4627-8a2c-d0878b31cc22
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ 12a3d153-4119-4fbd-af67-67e7227c8d97
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ 769b84d2-d3cb-4036-a979-37019b414eb6
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
## Tabu Search
"""

# ╔═╡ fab77d3b-203d-4033-bdd0-6ee1ff3cc572
function mhTS(tour::Array{Int,1}, dist::Array{Float64,2}, nbh::Array{Function,1})
    n, _ = size(dist)
    ctour = copy(tour)
    besttour = copy(tour)

    tabulist::Array{Array{Int,1},1} = []
    tabutenure = 7
    maxiter = 10 * n

    busqueda0 = [tspDist(dist, ctour)] # actual ztour
    busqueda1 = copy(busqueda0)        # mejor  ztour
    busquedafn = ["base"]              # actual fn

    for _ in 1:maxiter
        nbh_tour = [nbh_i(ctour, dist, tabulist) for nbh_i in nbh]
        nbh_ztour = []
        for (t, mv1, fn) in nbh_tour
            if t === nothing
                continue
            end
            zt = tspDist(dist, t)
            push!(nbh_ztour, (zt, t, mv1, fn))
        end
        if length(nbh_ztour) == 0
            break
        end

        sort!(nbh_ztour)
        newtourz, newtour, mov1, fn = nbh_ztour[1]

        ctour = newtour
        if newtourz < busqueda1[end]
            besttour = newtour
        end

        for mvi in mov1
            push!(tabulist, mvi)
        end
        while length(tabulist) > tabutenure
            popfirst!(tabulist)
        end

        push!(busqueda0, newtourz)
        if newtourz < busqueda1[end]
            push!(busqueda1, newtourz)
        else
            push!(busqueda1, busqueda1[end])
        end
        push!(busquedafn, fn)

        @show newtourz, fn
    end

    plot(busqueda0)
    plot!(busqueda1)
    savefig(joinpath(@__DIR__, "tsp-plot-TSz.png"))

    return besttour
end

# ╔═╡ cc590562-b6d7-422a-bdfd-267accb8d14a
# nbh::Array{Function,1} = [tspInsTS]
# nbh::Array{Function,1} = [tspSwapTS]
# nbh::Array{Function,1} = [tsp2OptTS]
# nbh::Array{Function,1} = [tspInsTS, tspSwapTS, tsp2OptTS]
# nbh::Array{Function,1} = [tsp2OptTS, tspSwapTS, tspInsTS]
nbh::Array{Function,1} = [tspSwapTS, tspInsTS, tsp2OptTS]

# ╔═╡ fe14fe38-c519-4103-bf59-207cfd2eb402
# tourMH = mhTS(tour, dist, nbh)
tourMH = mhTS(tour, dist, nbh)

# ╔═╡ f5050494-1975-4558-9a1f-2a44b388b3dd
tspPlot(X, Y, tourMH)

# ╔═╡ 0c85af3d-16bc-493c-a4e2-43af9a61026c
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-TSz.png"))

# ╔═╡ d14ed81d-f3e6-4067-9f19-faa93d6047a5
begin
    # using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ Cell order:
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═5329815c-850b-4d8b-9846-a4ffdf02b448
# ╠═879c42ee-f980-42e9-bb81-2f27672aecbf
# ╠═03fcd07b-99d5-4e2d-bc7f-9e00909139c6
# ╠═23121b1a-e4b1-4b47-8a2a-05fb11640818
# ╠═e03fba74-b7f1-4203-a8a2-faa795a514b8
# ╠═da135245-37a5-40a1-b29b-2ba403f3af47
# ╠═28e27c70-5d42-4096-a40d-ba64848865e4
# ╠═3817cb46-09bb-41ff-8dd1-d345bc8f5136
# ╠═a30695c5-2057-4627-8a2c-d0878b31cc22
# ╠═12a3d153-4119-4fbd-af67-67e7227c8d97
# ╠═769b84d2-d3cb-4036-a979-37019b414eb6
# ╠═d049dd39-4736-4ff0-b2b0-981c61ee3296
# ╠═fab77d3b-203d-4033-bdd0-6ee1ff3cc572
# ╠═cc590562-b6d7-422a-bdfd-267accb8d14a
# ╠═fe14fe38-c519-4103-bf59-207cfd2eb402
# ╠═f5050494-1975-4558-9a1f-2a44b388b3dd
# ╠═0c85af3d-16bc-493c-a4e2-43af9a61026c
# ╠═d14ed81d-f3e6-4067-9f19-faa93d6047a5

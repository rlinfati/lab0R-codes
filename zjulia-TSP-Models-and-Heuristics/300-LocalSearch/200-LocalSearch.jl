### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 413d6f34-d13a-47dd-b711-98d81f50f55e
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

# ╔═╡ 51bfede6-1436-11ec-0cc6-bfff507e0e99
using Plots

# ╔═╡ 03fcd07b-99d5-4e2d-bc7f-9e00909139c6
using Random

# ╔═╡ 07e5dfcb-ecb1-4310-87cd-5647cd3929c2
using PlutoUI

# ╔═╡ 263d3dc2-23ef-43a0-a38f-693abc3ff6fe
include(joinpath(@__DIR__, "../200-Heuristic/100-function.jl"))

# ╔═╡ a732af47-4eab-4d7c-9727-dc03a9143e29
include(joinpath(@__DIR__, "../200-Heuristic/200-initial.jl"))

# ╔═╡ fbec0e7d-2d67-4e8c-ad41-b42a5dd2e52e
include(joinpath(@__DIR__, "../300-LocalSearch/100-neighborhood.jl"))

# ╔═╡ 0e768f80-b14a-4095-9bd0-2a7523715e83
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 4a6a7189-6f6e-4211-b2fc-5e5b48aec3f9
n = 7

# ╔═╡ b03019d3-0270-45d8-ab5e-fcfb8d20dbc6
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ 3f648c67-e281-406e-8953-28c6f138280c
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ 1f47d25b-69b6-4c18-a1c1-5ce8ef922611
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

# ╔═╡ 4a1c0d9a-f99d-4402-ac9f-c341d04164de
tspDist(dist, tour)

# ╔═╡ 8d24ee93-ab92-49f6-bcfe-10a491f0a9d9
tspPlot(X, Y, tour)

# ╔═╡ d049dd39-4736-4ff0-b2b0-981c61ee3296
md"""
## Busqueda Local - Un Vecindario
"""

# ╔═╡ fab77d3b-203d-4033-bdd0-6ee1ff3cc572
function mhLocalSearch1(tour::Array{Int,1}, dist::Array{Float64,2}, nbh::Function; sf::Bool = false)
    ctour = copy(tour)
    busqueda0 = [tspDist(dist, ctour)] # actual ztour

    while true
        newtour = nbh(ctour, dist, sf)
        if newtour === nothing
            break
        end

        newtourz = tspDist(dist, newtour)
        @show newtourz, newtour

        if newtourz < busqueda0[end]
            push!(busqueda0, newtourz)
            ctour = newtour
        else
            break
        end
    end

    savefig(plot(busqueda0), joinpath(@__DIR__, "tsp-plot-ls1.png"))
    return ctour
end

# ╔═╡ 3a6639ba-2c69-411d-aedf-a79f37af0bbf
# nbh1::Function = tspInsLS
# nbh1::Function = tspSwapLS
# nbh1::Function = tsp2OptLS
nbh1::Function = tsp2OptLS

# ╔═╡ 624930ac-db0b-4cad-9195-b7617287ecd0
# tourLS1 = mhLocalSearch1(tour, dist, nbh1)
# tourLS1 = mhLocalSearch1(tour, dist, nbh1; sf = true)
# tourLS1 = mhLocalSearch1(tour, dist, nbh1; sf = false)
tourLS1 = mhLocalSearch1(tour, dist, nbh1, sf = true)

# ╔═╡ d8b2fc58-194a-46dd-bd0e-c85603719971
tspPlot(X, Y, tourLS1)

# ╔═╡ a2340e34-f4d7-41ee-b923-b6b2074fb3e5
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-ls1.png"))

# ╔═╡ 12216cdb-e83b-4e19-825c-c029a6781319
md"""
## Busqueda Local - Varios Vecindarios
"""

# ╔═╡ 35265c1c-11a2-457c-a005-ff96cc31eb89
function mhLocalSearch(tour::Array{Int,1}, dist::Array{Float64,2}, nbh::Array{Function,1}; sf::Bool = false)
    ctour = copy(tour)
    busqueda0 = [tspDist(dist, ctour)] # actual ztour
    busquedafn = ["base"] # actual fn

    while true
        nbh_tour = [(nbh_i(ctour, dist, sf), string(nbh_i)) for nbh_i in nbh]
        nbh_ztour = []
        for (t, fn) in nbh_tour
            if t === nothing
                continue
            end
            dzt = tspDist(dist, t) - busqueda0[end]
            if dzt > 0.0
                continue
            end
            push!(nbh_ztour, (dzt, t, fn))
        end

        if length(nbh_ztour) == 0
            break
        end

        if sf == false
            sort!(nbh_ztour)
        end

        ctour = nbh_ztour[1][2]
        push!(busqueda0, tspDist(dist, ctour))
        push!(busquedafn, nbh_ztour[1][3])

        _z = busqueda0[end]
        _f = busquedafn[end]

        @show _z, _f, ctour
    end

    savefig(plot(busqueda0), joinpath(@__DIR__, "tsp-plot-lss.png"))
    return ctour
end

# ╔═╡ fe14fe38-c519-4103-bf59-207cfd2eb402
# nbh::Array{Function,1} = [tspInsLS]
# nbh::Array{Function,1} = [tspSwapLS]
# nbh::Array{Function,1} = [tsp2OptLS]
# nbh::Array{Function,1} = [tspInsLS, tspSwapLS, tsp2OptLS]
# nbh::Array{Function,1} = [tsp2OptLS, tspSwapLS, tspInsLS]
nbh::Array{Function,1} = [tspSwapLS, tspInsLS]

# ╔═╡ 0dc3bc33-04c7-4905-8245-78e4412a47fd
# tourLSS = mhLocalSearch(tour, dist, nbh)
# tourLSS = mhLocalSearch(tour, dist, nbh; sf = true)
# tourLSS = mhLocalSearch(tour, dist, nbh; sf = false)
tourLSS = mhLocalSearch(tour, dist, nbh, sf = true)

# ╔═╡ f5050494-1975-4558-9a1f-2a44b388b3dd
tspPlot(X, Y, tourLSS)

# ╔═╡ 1faf1c9d-9717-45e5-9950-91400383c71b
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-lss.png"))

# ╔═╡ 35fc2435-c3aa-4c9b-ab60-4418c1b6a91a
begin
    # using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ Cell order:
# ╠═413d6f34-d13a-47dd-b711-98d81f50f55e
# ╠═0e768f80-b14a-4095-9bd0-2a7523715e83
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═03fcd07b-99d5-4e2d-bc7f-9e00909139c6
# ╠═07e5dfcb-ecb1-4310-87cd-5647cd3929c2
# ╠═263d3dc2-23ef-43a0-a38f-693abc3ff6fe
# ╠═a732af47-4eab-4d7c-9727-dc03a9143e29
# ╠═fbec0e7d-2d67-4e8c-ad41-b42a5dd2e52e
# ╠═4a6a7189-6f6e-4211-b2fc-5e5b48aec3f9
# ╠═b03019d3-0270-45d8-ab5e-fcfb8d20dbc6
# ╠═3f648c67-e281-406e-8953-28c6f138280c
# ╠═1f47d25b-69b6-4c18-a1c1-5ce8ef922611
# ╠═4a1c0d9a-f99d-4402-ac9f-c341d04164de
# ╠═8d24ee93-ab92-49f6-bcfe-10a491f0a9d9
# ╠═d049dd39-4736-4ff0-b2b0-981c61ee3296
# ╠═fab77d3b-203d-4033-bdd0-6ee1ff3cc572
# ╠═3a6639ba-2c69-411d-aedf-a79f37af0bbf
# ╠═624930ac-db0b-4cad-9195-b7617287ecd0
# ╠═d8b2fc58-194a-46dd-bd0e-c85603719971
# ╠═a2340e34-f4d7-41ee-b923-b6b2074fb3e5
# ╠═12216cdb-e83b-4e19-825c-c029a6781319
# ╠═35265c1c-11a2-457c-a005-ff96cc31eb89
# ╠═fe14fe38-c519-4103-bf59-207cfd2eb402
# ╠═0dc3bc33-04c7-4905-8245-78e4412a47fd
# ╠═f5050494-1975-4558-9a1f-2a44b388b3dd
# ╠═1faf1c9d-9717-45e5-9950-91400383c71b
# ╠═35fc2435-c3aa-4c9b-ab60-4418c1b6a91a

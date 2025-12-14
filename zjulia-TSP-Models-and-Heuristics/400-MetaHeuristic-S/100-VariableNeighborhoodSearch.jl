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

# ╔═╡ 47d0355d-06ce-4d73-bfea-1f2451386f32
using PlutoUI

# ╔═╡ ac82a645-a0c4-42a7-af5c-272105943f5a
include(joinpath(@__DIR__, "../200-Heuristic/100-function.jl"))

# ╔═╡ bb165467-c247-47eb-b186-cf234c279da6
include(joinpath(@__DIR__, "../200-Heuristic/200-initial.jl"))

# ╔═╡ 3ebea908-ef12-4ad8-a166-dcd702a5e3e7
include(joinpath(@__DIR__, "../300-LocalSearch/100-neighborhood.jl"))

# ╔═╡ ad7e3229-7b78-4eba-9beb-dfa0e8dd3b63
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 99c5e9cb-fb58-4bd2-990f-f28a99e30174
n = 7

# ╔═╡ 78572a1e-aed7-4c53-98a0-1e06b0ec0d22
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ 524cf8c9-3f4c-4738-8f5e-32b198ad8ab8
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ cb1f9214-a06d-46ed-9edd-55dc32954ef2
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

# ╔═╡ 90497267-22a8-400b-9799-9b0dfb629f05
md"""
## Variable Neighborhood Search
"""

# ╔═╡ 7065b5f0-b6a5-49ad-9dc4-08bf704e50e6
function mhVNS(tour::Array{Int,1}, dist::Array{Float64,2}, nbh::Array{Function,1}; sf::Bool = false)
    ctour = copy(tour)
    busqueda0 = [tspDist(dist, ctour)] # actual ztour
    busquedafn = ["base"] # actual fn

    while true
        nomejora = false
        for nbh_i in nbh
            newtour = nbh_i(ctour, dist, sf)
            if nbh_i == nbh[end]
                nomejora = true
            end
            if newtour === nothing
                continue
            end

            zt = tspDist(dist, newtour)
            if zt < busqueda0[end]
                ctour = newtour
                push!(busqueda0, zt)
                push!(busquedafn, string(nbh_i))
                @show zt, nbh_i, ctour
                nomejora = false
                break
            end
        end
        if nomejora
            break
        end
    end

    savefig(plot(busqueda0), joinpath(@__DIR__, "tsp-plot-vns.png"))

    return ctour
end

# ╔═╡ f81a0464-968d-4b8e-9688-bfecdd086f80
# nbh::Array{Function,1} = [tspInsLS]
# nbh::Array{Function,1} = [tspSwapLS]
# nbh::Array{Function,1} = [tsp2OptLS]
# nbh::Array{Function,1} = [tspInsLS, tspSwapLS, tsp2OptLS]
# nbh::Array{Function,1} = [tsp2OptLS, tspSwapLS, tspInsLS]
nbhvns::Array{Function,1} = [tspSwapLS, tspInsLS, tsp2OptLS]

# ╔═╡ 98d6408b-2d04-42bd-aecb-2ad3f9ff742b
# tourMH = mhVNS(tour, dist, nbh)
# tourMH = mhVNS(tour, dist, nbh; sf = true)
# tourMH = mhVNS(tour, dist, nbh; sf = false)
tourMH = mhVNS(tour, dist, nbhvns; sf = true)

# ╔═╡ fd0cda02-86c6-47c8-aaf2-5ffa9d6f2f6b
tspPlot(X, Y, tourMH)

# ╔═╡ b626019a-5639-4c12-bd28-4b0da0ae2065
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-vns.png"))

# ╔═╡ 16d6e720-c9f8-4d75-bc58-5eb505c738c3
begin
    # using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ Cell order:
# ╠═413d6f34-d13a-47dd-b711-98d81f50f55e
# ╠═ad7e3229-7b78-4eba-9beb-dfa0e8dd3b63
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═03fcd07b-99d5-4e2d-bc7f-9e00909139c6
# ╠═47d0355d-06ce-4d73-bfea-1f2451386f32
# ╠═ac82a645-a0c4-42a7-af5c-272105943f5a
# ╠═bb165467-c247-47eb-b186-cf234c279da6
# ╠═3ebea908-ef12-4ad8-a166-dcd702a5e3e7
# ╠═99c5e9cb-fb58-4bd2-990f-f28a99e30174
# ╠═78572a1e-aed7-4c53-98a0-1e06b0ec0d22
# ╠═524cf8c9-3f4c-4738-8f5e-32b198ad8ab8
# ╠═cb1f9214-a06d-46ed-9edd-55dc32954ef2
# ╠═4a1c0d9a-f99d-4402-ac9f-c341d04164de
# ╠═8d24ee93-ab92-49f6-bcfe-10a491f0a9d9
# ╠═90497267-22a8-400b-9799-9b0dfb629f05
# ╠═7065b5f0-b6a5-49ad-9dc4-08bf704e50e6
# ╠═f81a0464-968d-4b8e-9688-bfecdd086f80
# ╠═98d6408b-2d04-42bd-aecb-2ad3f9ff742b
# ╠═fd0cda02-86c6-47c8-aaf2-5ffa9d6f2f6b
# ╠═b626019a-5639-4c12-bd28-4b0da0ae2065
# ╠═16d6e720-c9f8-4d75-bc58-5eb505c738c3

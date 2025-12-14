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

# ╔═╡ 9b8872d8-9b2d-4a1f-82e2-2e536d8912ba
using Plots

# ╔═╡ d7286808-faad-4652-b2e2-27b77bf63e61
using Random

# ╔═╡ 03fcd07b-99d5-4e2d-bc7f-9e00909139c6
using Statistics

# ╔═╡ b48c6ba8-6302-4384-ac5f-534b85044a90
using PlutoUI

# ╔═╡ 133534ec-d6a9-43f5-a503-bf6126c8a865
include(joinpath(@__DIR__, "../200-Heuristic/100-function.jl"))

# ╔═╡ 5479487a-0fc7-4d26-94d8-8406b48aebca
include(joinpath(@__DIR__, "../200-Heuristic/200-initial.jl"))

# ╔═╡ ea29fcd9-ac2d-4b13-8809-e8925871184b
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 6ef09464-2e99-43c3-99e9-af732d04db80
md"""
## Ant Colony Optimization
"""

# ╔═╡ 8c2906e8-9656-44d2-89c4-0a5c2061440f
function mhACO(dist::Array{Float64,2}; s::Int = 0)
    n, _ = size(dist)
    if s == 0
        s = rand(1:n)
    end
    @assert 1 <= s <= n

    alpha = 1.0
    betha = 1.0
    rho = 0.20
    qfer = 100.0

    nHormigas = 2 * n
    maxiter = 10 * n

    tau::Array{Float64,2} = ones(n, n)
    eta::Array{Float64,2} = mean(dist) ./ dist
    for i in 1:n
        eta[i, i] = eps()
    end

    busqueda0 = [] # promedio ztour
    busqueda1 = [] # mejor    ztour
    mejortour = [0]

    for _ in 1:maxiter
        tours = []
        for _ in 1:nHormigas
            tour = Int[s]

            novisited = ones(n)
            novisited[s] = 0.0

            while sum(novisited) > 0.0
                pxy = tau[tour[end], :] .^ alpha .* eta[tour[end], :] .^ betha .* novisited
                pxy = pxy ./ sum(pxy)

                next = 0
                cspxy = cumsum(pxy)
                r = rand()

                for i in 1:n
                    if novisited[i] ≈ 0.0
                        continue
                    end
                    if r < cspxy[i]
                        next = i
                        break
                    end
                end
                @assert next != 0

                push!(tour, next)
                novisited[next] = 0.0
            end
            push!(tour, tour[1])
            push!(tours, tour)
        end

        zt = [tspDist(dist, i) for i in tours]
        push!(busqueda0, mean(zt))
        push!(busqueda1, minimum(zt))

        bestID = argmin(zt)
        if zt[bestID] < tspDist(dist, mejortour)
            mejortour = tours[bestID]
        end

        tau *= 1 - rho
        for h in 1:nHormigas
            up = qfer / tspDist(dist, tours[h])
            for p in 1:n
                i = tours[h][p]
                j = tours[h][p+1]
                tau[i, j] += up
            end
        end
    end

    plot(busqueda0)
    plot!(busqueda1)
    savefig(joinpath(@__DIR__, "tsp-plot-ACOz.png"))

    return mejortour
end

# ╔═╡ efaaa6d9-e90a-42fd-8d82-560707d127fc
n = 7

# ╔═╡ 2e029e77-7bbe-4b48-b0b2-ffec11528aa7
# X, Y = tspXYRand(n)
# X, Y = tspXYCluster(n)
X, Y = tspXYRand(n)

# ╔═╡ c93e4f96-8169-4db2-a988-21fac0c1fd07
# dist = distEuclidean(X, Y)
# dist = distManhattan(X, Y)	
# dist = distMaximum(X, Y)
# dist = distGeographical(X, Y)
dist = distEuclidean(X, Y)

# ╔═╡ e25019e0-3c89-4dc9-8adc-96900affe702
tourMH = mhACO(dist)

# ╔═╡ f5050494-1975-4558-9a1f-2a44b388b3dd
tspPlot(X, Y, tourMH)

# ╔═╡ 1e099701-6fd7-42d3-a76e-0077c4c8a174
PlutoUI.LocalResource(joinpath(@__DIR__, "tsp-plot-ACOz.png"))

# ╔═╡ d88138ac-178d-457e-9bb8-ff141e8d7caf
begin
    # using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ Cell order:
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═ea29fcd9-ac2d-4b13-8809-e8925871184b
# ╠═9b8872d8-9b2d-4a1f-82e2-2e536d8912ba
# ╠═d7286808-faad-4652-b2e2-27b77bf63e61
# ╠═03fcd07b-99d5-4e2d-bc7f-9e00909139c6
# ╠═b48c6ba8-6302-4384-ac5f-534b85044a90
# ╠═133534ec-d6a9-43f5-a503-bf6126c8a865
# ╠═5479487a-0fc7-4d26-94d8-8406b48aebca
# ╠═6ef09464-2e99-43c3-99e9-af732d04db80
# ╠═8c2906e8-9656-44d2-89c4-0a5c2061440f
# ╠═efaaa6d9-e90a-42fd-8d82-560707d127fc
# ╠═2e029e77-7bbe-4b48-b0b2-ffec11528aa7
# ╠═c93e4f96-8169-4db2-a988-21fac0c1fd07
# ╠═e25019e0-3c89-4dc9-8adc-96900affe702
# ╠═f5050494-1975-4558-9a1f-2a44b388b3dd
# ╠═1e099701-6fd7-42d3-a76e-0077c4c8a174
# ╠═d88138ac-178d-457e-9bb8-ff141e8d7caf

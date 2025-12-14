### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 36184f24-de30-11ec-3a6f-6962a488690c
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    # Pkg.add([])
    # Pkg.status()
end

# ╔═╡ 1acdf778-7b93-493b-b6a6-5230bf9cda52
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ 58d22a16-eaba-49d6-98f9-a39ae8b15052
md"""
## Generación de Instancia 
"""

# ╔═╡ 478ba342-7945-4a58-b31f-e24cfe3286e7
function tspXYRand(n::Int; myseed::Int = 1234)
    rng = Random.MersenneTwister(myseed) # 1.6 compat
    X = rand(rng, n) * 1_000.0
    Y = rand(rng, n) * 1_000.0
    return X, Y
end

# ╔═╡ 921250b2-b216-49be-9ca5-47a4a27e6e52
function tspXYCluster(n::Int; myseed::Int = 1234)
    rng = Random.MersenneTwister(myseed) # 1.6 compat

    cs = rand(rng, 3:8)
    X = rand(rng, n) * 1_000.0
    Y = rand(rng, n) * 1_000.0

    for i in (cs+1):n
        while true
            dps = sum(exp(-sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) / 40.0) for j in 1:cs)

            if rand(rng) < dps
                break
            end

            X[i] = rand(rng) * 1_000.0
            Y[i] = rand(rng) * 1_000.0
        end
    end
    return X, Y
end

# ╔═╡ c21fb05d-149b-4fd5-993b-43830831141b
md"""
## Calculo de matriz de costos
"""

# ╔═╡ 78d82c87-cd42-49df-9254-d6b1c925db63
function distManhattan(X::Array{Float64,1}, Y::Array{Float64,1})
    n = length(X)
    @assert n == length(Y)
    dist = [abs(X[i] - X[j]) + abs(Y[i] - Y[j]) for i in 1:n, j in 1:n]
    return dist
end

# ╔═╡ 1db7e55f-2aee-4113-96d3-b20d2921792f
function distEuclidean(X::Array{Float64,1}, Y::Array{Float64,1})
    n = length(X)
    @assert n == length(Y)
    dist = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return dist
end

# ╔═╡ 48686b17-e2a9-43e8-8061-a4fc780f4656
function distMaximum(X::Array{Float64,1}, Y::Array{Float64,1})
    n = length(X)
    @assert n == length(Y)
    dist = [max(abs(X[i] - X[j]), abs(Y[i] - Y[j])) for i in 1:n, j in 1:n]
    return dist
end

# ╔═╡ a825a088-4d1f-4df3-8fdc-405d4c40816d
function distGeographical(X::Array{Float64,1}, Y::Array{Float64,1})
    n = length(X)
    @assert n == length(Y)
    # xy = DDD.MM format
    # dist = GEO from TSPLIB
    lat = floor.(X)
    lat += (X - lat) * 5.0 / 3.0
    lat *= π / 180.0

    lon = floor.(Y)
    lon += (Y - lon) * 5.0 / 3.0
    lon *= π / 180.0

    q1 = [cos(lon[i] - lon[j]) for i in 1:n, j in 1:n]
    q2 = [cos(lat[i] - lat[j]) for i in 1:n, j in 1:n]
    q3 = [cos(lat[i] + lat[j]) for i in 1:n, j in 1:n]

    rrr = 6378.388

    dist = 0.5 * ((1.0 .+ q1) .* q2 - (1.0 .- q1) .* q3)
    dist = rrr * acos.(dist)

    return dist
end

# ╔═╡ cfe0d607-d517-4bc6-9977-3138b70830a9
md"""
## Calculo de costos totales
"""

# ╔═╡ 7fd89ab8-f587-4c51-95ad-9d4f364aa838
function tspDist(dist::Array{Float64,2}, tour::Array{Int,1})
    n, n2 = size(dist)
    @assert n == n2
    @assert n >= maximum(tour)

    f = tour[1] == tour[end] && sort(tour[1:(end-1)]) == 1:n
    if f == false
        return +Inf
    end

    totalDist = sum(dist[tour[i], tour[i+1]] for i in 1:n)
    return totalDist
end

# ╔═╡ 1d0b4449-2c52-4c0a-bd08-b2a46b729d4d
md"""
### Solución gráfica
"""

# ╔═╡ ed400f92-4184-43b7-ac19-c60a527a8e92
function tspPlot(X::Array{Float64,1}, Y::Array{Float64,1}, tour::Array{Int,1})
    n = length(X)
    @assert n == length(Y)
    @assert length(tour) == n + 1

    p = plot(legend = false)
    scatter!(X, Y, color = :blue)
    for i in 1:n
        annotate!(X[i], Y[i], text("$i", :top))
    end

    plot!(X[tour], Y[tour], color = :black)

    return p
end

# ╔═╡ d3579179-182e-4ded-b8a9-4d9cc60aae33
md"""
## Solucion Inicial Random
"""

# ╔═╡ eed1099c-dc74-4301-9a14-3fc2b20d386f
function tspRND(dist::Array{Float64,2})
    n, _ = size(dist)

    tour = randperm(n)
    push!(tour, tour[1])

    return tour
end

# ╔═╡ f8d0e342-05c7-42d1-b2ea-188131f21d2a
nothing

# ╔═╡ Cell order:
# ╠═36184f24-de30-11ec-3a6f-6962a488690c
# ╠═1acdf778-7b93-493b-b6a6-5230bf9cda52
# ╠═58d22a16-eaba-49d6-98f9-a39ae8b15052
# ╠═478ba342-7945-4a58-b31f-e24cfe3286e7
# ╠═921250b2-b216-49be-9ca5-47a4a27e6e52
# ╠═c21fb05d-149b-4fd5-993b-43830831141b
# ╠═78d82c87-cd42-49df-9254-d6b1c925db63
# ╠═1db7e55f-2aee-4113-96d3-b20d2921792f
# ╠═48686b17-e2a9-43e8-8061-a4fc780f4656
# ╠═a825a088-4d1f-4df3-8fdc-405d4c40816d
# ╠═cfe0d607-d517-4bc6-9977-3138b70830a9
# ╠═7fd89ab8-f587-4c51-95ad-9d4f364aa838
# ╠═1d0b4449-2c52-4c0a-bd08-b2a46b729d4d
# ╠═ed400f92-4184-43b7-ac19-c60a527a8e92
# ╠═d3579179-182e-4ded-b8a9-4d9cc60aae33
# ╠═eed1099c-dc74-4301-9a14-3fc2b20d386f
# ╠═f8d0e342-05c7-42d1-b2ea-188131f21d2a

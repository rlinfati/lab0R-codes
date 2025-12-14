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
## Solucion Inicial
"""

# ╔═╡ 921250b2-b216-49be-9ca5-47a4a27e6e52
function tspNN(dist::Array{Float64,2}; s::Int = 0)
    n, n2 = size(dist)
    @assert n == n2

    if s == 0
        s = rand(1:n)
    end
    @assert 1 <= s <= n

    tour = Int[s]
    aVisitar = collect(1:n)
    deleteat!(aVisitar, s)

    while !isempty(aVisitar)
        proxID = argmin(dist[tour[end], aVisitar])
        push!(tour, aVisitar[proxID])
        deleteat!(aVisitar, proxID)
    end
    push!(tour, s)

    return tour
end

# ╔═╡ 7205e117-8f48-4958-b3b7-489d01bdcaf0
function tspCI(dist::Array{Float64,2}; s::Int = 0)
    n, n2 = size(dist)
    @assert n == n2

    if s == 0
        s = rand(1:n)
    end
    @assert 1 <= s <= n

    tour = Int[]
    push!(tour, s)
    push!(tour, s)

    aVisitar = collect(1:n)
    deleteat!(aVisitar, s)

    while !isempty(aVisitar)
        bestJ = 0
        bestK = 0

        f(idx, k) = dist[tour[idx], k] + dist[k, tour[idx+1]] - dist[tour[idx], tour[idx+1]]
        m = [f(idx, k) for idx in 1:(length(tour)-1), k in aVisitar]

        x = argmin(m)
        bestJ = x[1] + 1
        bestK = aVisitar[x[2]]

        insert!(tour, bestJ, bestK)
        setdiff!(aVisitar, bestK)
    end

    return tour
end

# ╔═╡ 3185c1e6-49a0-4033-bcdd-dcd025254f75
function tspCWs(dist::Array{Float64,2}; lamda::Float64 = 1.0, s::Int = 0)
    n, n2 = size(dist)
    @assert n == n2

    if s == 0
        s = rand(1:n)
    end
    @assert 1 <= s <= n

    ahorro = [(dist[i, s] + dist[s, j] - lamda * dist[i, j], i, j) for i in 1:n for j in (i+1):n]
    sort!(ahorro, rev = true)

    visitado = zeros(Int, n)

    tour = [ahorro[1][2]; ahorro[1][3]]

    visitado[s] = 1
    visitado[tour[1]] = 1
    visitado[tour[2]] = 1
    visitado[s] = 1

    while sum(visitado) < n
        for (_, i, j) in ahorro
            if visitado[i] + visitado[j] != 1
                continue
            end

            if tour[1] == i
                pushfirst!(tour, j)
                visitado[j] = 1
            elseif tour[1] == j
                pushfirst!(tour, i)
                visitado[i] = 1
            elseif tour[end] == i
                push!(tour, j)
                visitado[j] = 1
            elseif tour[end] == j
                push!(tour, i)
                visitado[i] = 1
            else
                continue
            end

            break
        end
    end

    pushfirst!(tour, s)
    push!(tour, s)

    return tour
end

# ╔═╡ 607f06e6-bebb-4a95-815d-66153732d22d
function tspCWp(dist::Array{Float64,2}; lamda::Float64 = 1.0, s::Int = 0)
    n, n2 = size(dist)
    @assert n == n2

    if s == 0
        s = rand(1:n)
    end
    @assert 1 <= s <= n

    ahorro = [(dist[i, s] + dist[s, j] - lamda * dist[i, j], i, j) for i in 1:n for j in (i+1):n]
    sort!(ahorro, rev = true)

    visitado = zeros(n)
    paths = []

    while sum(visitado) < n || length(paths) > 1
        for (_, i, j) in ahorro
            if visitado[i] + visitado[j] == 0
                push!(paths, [i, j])
                visitado[i] = visitado[j] = 1
            elseif visitado[i] + visitado[j] == 1
                for p in paths
                    if p[1] == i
                        pushfirst!(p, j)
                        visitado[j] = 1
                    elseif p[1] == j
                        pushfirst!(p, i)
                        visitado[i] = 1
                    elseif p[end] == i
                        push!(p, j)
                        visitado[j] = 1
                    elseif p[end] == j
                        push!(p, i)
                        visitado[i] = 1
                    end
                end
            elseif visitado[i] + visitado[j] == 2
                v1 = findfirst(p -> p[1] == i, paths)
                v1s = true
                if v1 === nothing
                    v1 = findfirst(p -> p[end] == i, paths)
                    v1s = false
                end
                if v1 === nothing
                    continue
                end

                v2 = findfirst(p -> p[1] == j, paths)
                v2s = true
                if v2 === nothing
                    v2 = findfirst(p -> p[end] == j, paths)
                    v2s = false
                end
                if v2 === nothing
                    continue
                end

                if v1 == v2
                    continue
                end

                vn = Int[]
                if v1s == true && v2s == true
                    vn = [reverse(paths[v1]); paths[v2]]
                elseif v1s == true && v2s == false
                    vn = [paths[v2]; paths[v1]]
                elseif v1s == false && v2s == true
                    vn = [paths[v1]; paths[v2]]
                elseif v1s == false && v2s == false
                    vn = [paths[v1]; reverse(paths[v2])]
                end
                push!(paths, vn)
                deleteat!(paths, max(v1, v2))
                deleteat!(paths, min(v1, v2))
            end
        end
    end

    tour = paths[1]
    push!(tour, tour[1])

    return tour
end

# ╔═╡ 35c8a75e-fdd1-447a-a506-d446175d1674
nothing

# ╔═╡ Cell order:
# ╠═36184f24-de30-11ec-3a6f-6962a488690c
# ╠═1acdf778-7b93-493b-b6a6-5230bf9cda52
# ╠═58d22a16-eaba-49d6-98f9-a39ae8b15052
# ╠═921250b2-b216-49be-9ca5-47a4a27e6e52
# ╠═7205e117-8f48-4958-b3b7-489d01bdcaf0
# ╠═3185c1e6-49a0-4033-bcdd-dcd025254f75
# ╠═607f06e6-bebb-4a95-815d-66153732d22d
# ╠═35c8a75e-fdd1-447a-a506-d446175d1674

### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 413d6f34-d13a-47dd-b711-98d81f50f55e
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    # Pkg.add([])
    # Pkg.status()
end

# ╔═╡ 51bfede6-1436-11ec-0cc6-bfff507e0e99
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
"""

# ╔═╡ c7c9a12f-f3ab-4251-888d-d1ac9b87aa60
md"""
## Busquedas Locales - Swap
"""

# ╔═╡ e44f264e-e87d-4498-b89a-cdd8519eddd1
function tspSwapLS(tour::Array{Int,1}, dist::Array{Float64,2}, sf::Bool)
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = 0.0
    for i in 2:n-1, j in i+2:n
        dd = 0.0
        dd += dist[tour[i-1], tour[j]] + dist[tour[j], tour[i+1]]
        dd += dist[tour[j-1], tour[i]] + dist[tour[i], tour[j+1]]
        dd -= dist[tour[i-1], tour[i]] + dist[tour[i], tour[i+1]]
        dd -= dist[tour[j-1], tour[j]] + dist[tour[j], tour[j+1]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
            if sf
                break
            end
        end
    end

    if besti == 0 || bestj == 0
        return nothing
    end

    newtour = copy(tour)
    newtour[besti], newtour[bestj] = newtour[bestj], newtour[besti]
    return newtour
end

# ╔═╡ 88db4b1e-1dbd-48db-893a-ce230f65c690
function tspSwapTS(tour::Array{Int,1}, dist::Array{Float64,2}, tabulist::Array{Array{Int,1},1})
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = +Inf
    for i in 2:n-1, j in i+2:n
        if [tour[i-1], tour[j]] in tabulist
            continue
        end
        if [tour[j], tour[i+1]] in tabulist
            continue
        end
        if [tour[j-1], tour[i]] in tabulist
            continue
        end
        if [tour[i], tour[j+1]] in tabulist
            continue
        end

        dd = 0.0
        dd += dist[tour[i-1], tour[j]] + dist[tour[j], tour[i+1]]
        dd += dist[tour[j-1], tour[i]] + dist[tour[i], tour[j+1]]
        dd -= dist[tour[i-1], tour[i]] + dist[tour[i], tour[i+1]]
        dd -= dist[tour[j-1], tour[j]] + dist[tour[j], tour[j+1]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
        end
    end

    if besti == 0 || bestj == 0
        return (nothing, nothing, nothing)
    end

    newtour = copy(tour)
    newtour[besti], newtour[bestj] = newtour[bestj], newtour[besti]

    mov1 = [tour[besti-1], tour[besti]]
    mov2 = [tour[besti], tour[besti+1]]
    mov3 = [tour[bestj-1], tour[bestj]]
    mov4 = [tour[bestj], tour[bestj+1]]

    return (newtour, [mov1, mov2, mov3, mov4], "tspSwapTS")
end

# ╔═╡ 4c5a0402-8ed5-45e9-8cbe-7ce5fb094655
function tspSwapSA(tour::Array{Int,1}, dist::Array{Float64,2})
    n, _ = size(dist)

    besti = rand(2:n-1)
    bestj = rand(besti+1:n)

    newtour = copy(tour)
    newtour[besti], newtour[bestj] = newtour[bestj], newtour[besti]
    return newtour
end

# ╔═╡ 6dadce5c-6f35-441a-a415-1b07c1e5b5ca
md"""
## Busquedas Locales - Insertion
"""

# ╔═╡ 8572781c-d545-4229-ab20-1fbd84678c7e
function tspInsLS(tour::Array{Int,1}, dist::Array{Float64,2}, sf::Bool)
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = 0.0
    for i in 2:n-1, j in i+2:n+1
        dd = 0.0
        dd -= dist[tour[i-1], tour[i]] + dist[tour[i], tour[i+1]]
        dd += dist[tour[i-1], tour[i+1]]
        dd -= dist[tour[j-1], tour[j]]
        dd += dist[tour[j-1], tour[i]] + dist[tour[i], tour[j]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
            if sf
                break
            end
        end
    end

    if besti == 0 || bestj == 0
        return nothing
    end

    newtour = copy(tour)
    insert!(newtour, bestj, newtour[besti])
    deleteat!(newtour, besti)

    return newtour
end

# ╔═╡ 45412885-0051-48d0-8920-357e9a1970c8
function tspInsTS(tour::Array{Int,1}, dist::Array{Float64,2}, tabulist::Array{Array{Int,1},1})
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = +Inf
    for i in 2:n-1, j in i+2:n+1
        if [tour[i-1], tour[i+1]] in tabulist
            continue
        end
        if [tour[j-1], tour[i]] in tabulist
            continue
        end
        if [tour[i], tour[j]] in tabulist
            continue
        end

        dd = 0.0
        dd -= dist[tour[i-1], tour[i]] + dist[tour[i], tour[i+1]]
        dd += dist[tour[i-1], tour[i+1]]
        dd -= dist[tour[j-1], tour[j]]
        dd += dist[tour[j-1], tour[i]] + dist[tour[i], tour[j]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
        end
    end

    if besti == 0 || bestj == 0
        return (nothing, nothing, nothing)
    end

    newtour = copy(tour)
    insert!(newtour, bestj, newtour[besti])
    deleteat!(newtour, besti)

    mov1 = [tour[besti-1], tour[besti]]
    mov2 = [tour[besti], tour[besti+1]]
    mov3 = [tour[bestj-1], tour[bestj]]

    return (newtour, [mov1, mov2, mov3], "tspInsTS")
end

# ╔═╡ b438deaf-fac0-425e-bf48-97823044ec11
function tspInsSA(tour::Array{Int,1}, dist::Array{Float64,2})
    n, _ = size(dist)

    besti = rand(2:n-1)
    bestj = rand(besti+2:n+1)

    newtour = copy(tour)
    insert!(newtour, bestj, newtour[besti])
    deleteat!(newtour, besti)

    return newtour
end

# ╔═╡ 16163f96-89ba-4bd1-869d-aa90c749024d
md"""
## Busquedas Locales - 2-Opt
"""

# ╔═╡ 64cdae57-1a43-4300-9c35-16e94104280f
function tsp2OptLS(tour::Array{Int,1}, dist::Array{Float64,2}, sf::Bool)
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = 0.0
    for i in 2:n-1, j in i+1:n
        dd = dist[tour[i-1], tour[j]] + dist[tour[i], tour[j+1]] - dist[tour[i-1], tour[i]] - dist[tour[j], tour[j+1]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
            if sf
                break
            end
        end
    end

    if besti == 0 || bestj == 0
        return nothing
    end

    return reverse(tour, besti, bestj)
end

# ╔═╡ 07dd4b37-9458-49ad-8638-77e919c730f3
function tsp2OptTS(tour::Array{Int,1}, dist::Array{Float64,2}, tabulist::Array{Array{Int,1},1})
    n, _ = size(dist)

    besti = 0
    bestj = 0
    bestd = +Inf
    for i in 2:n-1, j in i+1:n
        if [tour[i-1], tour[j]] in tabulist
            continue
        end
        if [tour[i], tour[j+1]] in tabulist
            continue
        end

        dd = dist[tour[i-1], tour[j]] + dist[tour[i], tour[j+1]] - dist[tour[i-1], tour[i]] - dist[tour[j], tour[j+1]]
        if dd < bestd - eps(Float16)
            besti = i
            bestj = j
            bestd = dd
        end
    end

    if besti == 0 || bestj == 0
        return (nothing, nothing, nothing)
    end

    newtour = reverse(tour, besti, bestj)

    mov1 = [tour[besti-1], tour[besti]]
    mov2 = [tour[bestj], tour[bestj+1]]

    return (newtour, [mov1, mov2], "tsp2OptTS")
end

# ╔═╡ 54457fa0-dcb4-42ee-98c5-67c4702a324e
function tsp2OptSA(tour::Array{Int,1}, dist::Array{Float64,2})
    n, _ = size(dist)

    besti = rand(2:n-1)
    bestj = rand(besti+1:n)

    return reverse(tour, besti, bestj)
end

# ╔═╡ f267c6c7-686f-49a3-9831-982e19b9e2f3
nothing

# ╔═╡ Cell order:
# ╠═413d6f34-d13a-47dd-b711-98d81f50f55e
# ╠═51bfede6-1436-11ec-0cc6-bfff507e0e99
# ╠═c7c9a12f-f3ab-4251-888d-d1ac9b87aa60
# ╠═e44f264e-e87d-4498-b89a-cdd8519eddd1
# ╠═88db4b1e-1dbd-48db-893a-ce230f65c690
# ╠═4c5a0402-8ed5-45e9-8cbe-7ce5fb094655
# ╠═6dadce5c-6f35-441a-a415-1b07c1e5b5ca
# ╠═8572781c-d545-4229-ab20-1fbd84678c7e
# ╠═45412885-0051-48d0-8920-357e9a1970c8
# ╠═b438deaf-fac0-425e-bf48-97823044ec11
# ╠═16163f96-89ba-4bd1-869d-aa90c749024d
# ╠═64cdae57-1a43-4300-9c35-16e94104280f
# ╠═07dd4b37-9458-49ad-8638-77e919c730f3
# ╠═54457fa0-dcb4-42ee-98c5-67c4702a324e
# ╠═f267c6c7-686f-49a3-9831-982e19b9e2f3

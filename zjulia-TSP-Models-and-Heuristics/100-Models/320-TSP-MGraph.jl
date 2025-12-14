### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ e7262522-ac65-11ec-0633-1d82420161db
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    Pkg.add(
        [
            Pkg.PackageSpec("JuMP")
            Pkg.PackageSpec("GLPK")
            Pkg.PackageSpec("Graphs")
            Pkg.PackageSpec("MetaGraphs")
            Pkg.PackageSpec("GraphPlot")
            Pkg.PackageSpec("Plots")
            Pkg.PackageSpec("PlutoUI")
        ],
    )
    Pkg.status()
end

# ╔═╡ e726254a-ac65-11ec-2f4c-1bead54e006a
using JuMP

# ╔═╡ a28468d5-b6fd-42d1-856f-4b9a8198a4fa
using GLPK

# ╔═╡ 70f14c5a-0809-426d-a8a0-7dee755da648
using Graphs

# ╔═╡ 94e5017d-0bf4-411c-b7f7-abdc8d5fc961
using MetaGraphs

# ╔═╡ 2f6763ba-366f-48d4-85c6-df7a12e57242
using GraphPlot

# ╔═╡ 04e5d936-545a-470e-93f3-6065b80a87fc
using Plots

# ╔═╡ bb4218e4-06f0-482b-9ab6-55fedd51d429
using Random

# ╔═╡ c4da83d1-4e75-45c6-b380-496f73fd48d9
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ cbbfc620-5416-48fa-8ffa-c2227153707e
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Gavish-Graves

Lectura recomendada:

```
Gavish, B., & Graves, S. (1981). Scheduling and routing in transportation and distribution systems: formulations and new relaxations.
```
"""

# ╔═╡ 12bdc286-b9da-4b79-b06d-d01106d11c16
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 66bea84e-7574-494c-abcf-eb8d092e778c
md"""
### Generación de Instancia
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
n = 7

# ╔═╡ f9cb6029-062b-436e-968b-9e8e4dd7b8f6
begin
    rng = Random.MersenneTwister(1234) # 1.6 compat
    X = rand(rng, n) * 1_000.0
    Y = rand(rng, n) * 1_000.0
    n < 11 ? [1:n X Y] : nothing
end

# ╔═╡ 72b74a0b-e001-445f-9c65-995c69854951
md"""
### Calculo de matriz de costos
"""

# ╔═╡ c2a2b9e2-170f-44eb-be42-ea3219836e91
c = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]

# ╔═╡ 7ce03388-585e-4b78-8dfe-52d2c8ba8602
md"""
### Grafo con MetaGraphs
"""

# ╔═╡ 4573d40c-957f-4834-829d-53943421b73d
G = MetaDiGraph(complete_digraph(n))

# ╔═╡ 7ccc5694-d49f-4ff1-b2b5-c1a8278d0433
for v in vertices(G)
    set_prop!(G, v, :name, 'a' + v - 1)
end

# ╔═╡ bd5d2775-75eb-4474-a906-fc553c08a0a6
for e in edges(G)
    set_prop!(G, e, :weight, c[e.src, e.dst])
end

# ╔═╡ 8906094b-5a56-4cb9-ac01-ae2da5cdea22
n < 11 ? [(v, get_prop(G, v, :name)) for v in vertices(G)] : nothing

# ╔═╡ 1dbce19b-79af-49c5-b366-740f1d350f80
n < 11 ? [(e.src, e.dst, get_prop(G, e, :weight)) for e in edges(G)] : nothing

# ╔═╡ 99280caa-ff78-44e9-93ee-d657ad578e9a
gplot(
    G,
    X,
    -Y,
    nodelabel = [get_prop(G, v, :name) for v in vertices(G)],
    edgelabel = [get_prop(G, e, :weight) for e in edges(G)],
)

# ╔═╡ 495ebd29-608e-484d-904a-d21c869a7791
md"""
### Modelo en JuMP
"""

# ╔═╡ 249acea1-9778-4e03-9a7e-5636054c3f06
begin
    Node = vertices(G)
    Nod2 = setdiff(Node, first(Node))
    Arcs = edges(G)
end

# ╔═╡ 3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
begin
    m = JuMP.Model()

    @variable(m, x[Arcs], Bin)
    @variable(m, f[Arcs] >= 0)

    @objective(m, Min, sum(get_prop(G, e, :weight) * x[e] for e in Arcs))

    @constraint(m, r1[i in Node], sum(x[e] for e in Arcs if e.src == i) == 1)
    @constraint(m, r2[j in Node], sum(x[e] for e in Arcs if e.dst == j) == 1)

    # SEC Gavish and Graves formulation
    @constraint(m, r3, sum(f[e] for e in Arcs if e.src == 1) == n)
    @constraint(m, r4, sum(f[e] for e in Arcs if e.dst == 1) == 1)
    @constraint(m, r5[e in Arcs], x[e] <= f[e])
    @constraint(m, r6[e in Arcs], f[e] <= n * x[e])
    @constraint(m, r7[k in Nod2], sum(f[e] for e in Arcs if e.dst == k) == 1 + sum(f[e] for e in Arcs if e.src == k))

    n < 11 ? JuMP.latex_formulation(m) : nothing
end

# ╔═╡ 2edd33de-ebd4-484c-a45b-3fe7c6beb597
md"""
### Parametros del Solver y Solución
"""

# ╔═╡ d1b599d9-da3a-4070-90dc-e63532951fd6
begin
    JuMP.set_optimizer(m, GLPK.Optimizer)
    JuMP.set_optimizer_attribute(m, "msg_lev", GLP_MSG_ALL)
    JuMP.set_optimizer_attribute(m, "tm_lim", 60 * 1000)
    JuMP.optimize!(m)
end

# ╔═╡ fc925cdf-ffd8-45ad-a7a9-4c11228fac02
JuMP.solution_summary(m)

# ╔═╡ 21cc515d-d381-445f-a23a-4bc75c81e38c
md"""
### Solución del Solver
"""

# ╔═╡ 2fd0b79d-4740-45aa-a63d-95d2fb41ad81
xval = JuMP.value.(x) .≈ 1.0

# ╔═╡ 24683b46-ecc8-487c-81ec-5dff1e28dd87
fval = round.(Int, JuMP.value.(f))

# ╔═╡ 8700444f-d4b6-49ce-8551-14c877793ffc
md"""
### Solución gráfica
"""

# ╔═╡ 8a08b116-152a-40d5-9c28-b1a0c3524aa6
t = [e for e in Arcs if xval[e]]

# ╔═╡ 773c8409-3e01-4944-8ba4-0442b8a2368d
sum([get_prop(G, e, :weight) for e in t])

# ╔═╡ c2c05477-7c33-480f-b88f-c173487a4a4c
begin
    p = plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in Node
        annotate!(X[i], Y[i], text(get_prop(G, i, :name), :top))
    end

    for e in t
        e = [e.src; e.dst]
        plot!(X[e], Y[e], color = :black)
    end

    p
end

# ╔═╡ d6bf3ba9-3fa4-4257-bab0-1f2605a6585e
md"""
### Calculo del vector tour
"""

# ╔═╡ c806912b-8bec-4d41-a7f3-0c7654a4ed53
begin
    tour = Int[1]
    while true
        push!(tour, filter(e -> e.src == tour[end], t)[1].dst)
        if tour[end] == 1
            break
        end
    end
    tour
end

# ╔═╡ Cell order:
# ╠═e7262522-ac65-11ec-0633-1d82420161db
# ╠═cbbfc620-5416-48fa-8ffa-c2227153707e
# ╠═12bdc286-b9da-4b79-b06d-d01106d11c16
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═70f14c5a-0809-426d-a8a0-7dee755da648
# ╠═94e5017d-0bf4-411c-b7f7-abdc8d5fc961
# ╠═2f6763ba-366f-48d4-85c6-df7a12e57242
# ╠═04e5d936-545a-470e-93f3-6065b80a87fc
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═66bea84e-7574-494c-abcf-eb8d092e778c
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═7ce03388-585e-4b78-8dfe-52d2c8ba8602
# ╠═4573d40c-957f-4834-829d-53943421b73d
# ╠═7ccc5694-d49f-4ff1-b2b5-c1a8278d0433
# ╠═bd5d2775-75eb-4474-a906-fc553c08a0a6
# ╠═8906094b-5a56-4cb9-ac01-ae2da5cdea22
# ╠═1dbce19b-79af-49c5-b366-740f1d350f80
# ╠═99280caa-ff78-44e9-93ee-d657ad578e9a
# ╠═495ebd29-608e-484d-904a-d21c869a7791
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═2edd33de-ebd4-484c-a45b-3fe7c6beb597
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╟─21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╠═24683b46-ecc8-487c-81ec-5dff1e28dd87
# ╠═8700444f-d4b6-49ce-8551-14c877793ffc
# ╠═8a08b116-152a-40d5-9c28-b1a0c3524aa6
# ╠═773c8409-3e01-4944-8ba4-0442b8a2368d
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═d6bf3ba9-3fa4-4257-bab0-1f2605a6585e
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═c4da83d1-4e75-45c6-b380-496f73fd48d9

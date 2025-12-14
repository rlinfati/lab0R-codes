### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ f32a7be9-48a7-4195-b313-d1ed1219ea31
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    Pkg.add([
        Pkg.PackageSpec("JuMP")
        Pkg.PackageSpec("GLPK")
        Pkg.PackageSpec("Plots")
        Pkg.PackageSpec("PlutoUI")
    ])
    Pkg.status()
end

# ╔═╡ e726254a-ac65-11ec-2f4c-1bead54e006a
using JuMP

# ╔═╡ a28468d5-b6fd-42d1-856f-4b9a8198a4fa
using GLPK

# ╔═╡ e6902e75-32f3-4c61-b564-24ce502f7025
using Plots

# ╔═╡ bb4218e4-06f0-482b-9ab6-55fedd51d429
using Random

# ╔═╡ 3d44932a-5575-4663-ae0a-0c8cee944397
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ 55c4763f-55cd-469a-8900-c51bdc2c108a
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Gavish-Graves

Lectura recomendada:

```
Gavish, B., & Graves, S. (1981). Scheduling and routing in transportation and distribution systems: formulations and new relaxations.
```
"""

# ╔═╡ 5435d1a7-d8d8-4460-9a33-0665cb5d7276
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ f11ae2b0-123b-4ebf-8e6f-12ddce90231c
md"""
### Generación de Instancia
"""

# ╔═╡ 6625de43-fc6f-4545-9793-beabdef03426
n = 7

# ╔═╡ f9cb6029-062b-436e-968b-9e8e4dd7b8f6
begin
    rng = Random.MersenneTwister(1234) # 1.6 compat
    X = rand(rng, n) * 1_000.0
    Y = rand(rng, n) * 1_000.0
    n < 11 ? [1:n X Y] : nothing
end

# ╔═╡ a27390b6-39e3-422e-8832-a808f0a1b812
md"""
### Calculo de matriz de costos
"""

# ╔═╡ 72b74a0b-e001-445f-9c65-995c69854951
begin
    c = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    n < 11 ? c : nothing
end

# ╔═╡ 3bdf9fa5-3776-48de-a196-73e169652017
md"""
### Modelo en JuMP
"""

# ╔═╡ 4a26798d-737d-45d9-92b0-9fecaba005f8
begin
    Node = collect(1:n)
    Nod2 = collect(2:n)
    Arcs = [(o = i, d = j) for i in 1:n, j in 1:n if i != j]
    Dist = [c[e[:o], e[:d]] for e in Arcs]
    nothing
end

# ╔═╡ 3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
begin
    m = JuMP.Model()

    @variable(m, x[Arcs], Bin)
    @variable(m, f[Arcs] >= 0)

    @objective(m, Min, Dist' * x)

    @constraint(m, r1[i in Node], sum(x[e] for e in Arcs if e[:o] == i) == 1)
    @constraint(m, r2[j in Node], sum(x[e] for e in Arcs if e[:d] == j) == 1)

    # SEC Gavish and Graves formulation
    @constraint(m, r3, sum(f[e] for e in Arcs if e[:o] == 1) == n)
    @constraint(m, r4, sum(f[e] for e in Arcs if e[:d] == 1) == 1)
    @constraint(m, r5[e in Arcs], x[e] <= f[e])
    @constraint(m, r6[e in Arcs], f[e] <= n * x[e])
    @constraint(m, r7[k in Nod2], sum(f[e] for e in Arcs if e[:d] == k) == 1 + sum(f[e] for e in Arcs if e[:o] == k))

    n < 11 ? JuMP.latex_formulation(m) : nothing
end

# ╔═╡ 48e6f211-27a6-42a2-8c4a-4352863da1fa
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

# ╔═╡ 9da48217-e40c-40f1-a776-c8378e393aab
md"""
### Solución gráfica
"""

# ╔═╡ 8a08b116-152a-40d5-9c28-b1a0c3524aa6
t = [e for e in Arcs if xval[e]]

# ╔═╡ 33c830a1-9fdc-430a-85c0-988723569738
sum(c[e[:o], e[:d]] for e in t)

# ╔═╡ c2c05477-7c33-480f-b88f-c173487a4a4c
begin
    p = plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in Node
        annotate!(X[i], Y[i], text("$i", :top))
    end

    for e in t
        e = [e[:o]; e[:d]]
        plot!(X[e], Y[e], color = :black)
    end

    p
end

# ╔═╡ d6bf3ba9-3fa4-4257-bab0-1f2605a6585e
md"""
### Calculo del vector tour
"""

# ╔═╡ c782e93f-d22d-4292-82f6-4c2123dda7f2
sort!(t)

# ╔═╡ c806912b-8bec-4d41-a7f3-0c7654a4ed53
begin
    tour = Int[1]
    while true
        push!(tour, filter(e -> e[:o] == tour[end], t)[1][:d])
        if tour[end] == 1
            break
        end
    end
    tour
end

# ╔═╡ Cell order:
# ╠═f32a7be9-48a7-4195-b313-d1ed1219ea31
# ╠═55c4763f-55cd-469a-8900-c51bdc2c108a
# ╠═5435d1a7-d8d8-4460-9a33-0665cb5d7276
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═f11ae2b0-123b-4ebf-8e6f-12ddce90231c
# ╠═6625de43-fc6f-4545-9793-beabdef03426
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═a27390b6-39e3-422e-8832-a808f0a1b812
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═3bdf9fa5-3776-48de-a196-73e169652017
# ╠═4a26798d-737d-45d9-92b0-9fecaba005f8
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═48e6f211-27a6-42a2-8c4a-4352863da1fa
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╠═21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╠═24683b46-ecc8-487c-81ec-5dff1e28dd87
# ╠═9da48217-e40c-40f1-a776-c8378e393aab
# ╠═8a08b116-152a-40d5-9c28-b1a0c3524aa6
# ╠═33c830a1-9fdc-430a-85c0-988723569738
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═d6bf3ba9-3fa4-4257-bab0-1f2605a6585e
# ╠═c782e93f-d22d-4292-82f6-4c2123dda7f2
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═3d44932a-5575-4663-ae0a-0c8cee944397

### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ e7262522-ac65-11ec-0633-1d82420161db
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

# ╔═╡ 8b45efd6-0fec-4a1f-a5a6-2fa8d95e9a99
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ c22a2d8c-885b-4b52-a411-4eaf6d29c059
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Gavish-Graves

Lectura recomendada:

```
Gavish, B., & Graves, S. (1981). Scheduling and routing in transportation and distribution systems: formulations and new relaxations.
```
"""

# ╔═╡ 791783f2-1f1f-44cc-bd7e-8425e0d18c2f
md"""
## Modelo Matemático

Parametro $n$: Cantidad de ciudades.

Parametro $c_{ij}$: Costo (distancia, tiempo, etc.) de recorrer desde la ciudad $i$ a la ciudad $j$. Generalmente distancia euclidiana.

Variable $x_{ij}$ : Variable binaria que toma el valor 1 si el vendedor recorrer desde la ciudad $i$ a la ciudad $j$. 0 en otro caso (eoc).

Variable $f_{ij}$ : Variable continua, si el arco es visitado ($x_{ij} = 1$) etiqueta las ciudades con en orden inverso de visita.

$Min \: z = \sum_{i=1}^{n} \sum_{\displaylines{j=1 \\ i \neq j}}^{n} c_{ij}x_{ij}$

$\sum_{\displaylines{i=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall j = 1..n$

$\sum_{\displaylines{j=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall i = 1..n$

$\sum_{j=1}^{n} f_{1j} = n$

$\sum_{i=1}^{n} f_{i1} = 1$

$x_{ij} \leq f_{ij} \leq n x_{ij} \qquad i,j = 1..n$

$\sum_{i=1}^{n} f_{ik} = 1 + \sum_{j=1}^{n} f_{kj} \qquad k = 2..n$

$x_{ij} \in \{0,1\} \qquad \forall i,j=1..n; i \neq j$

$f_{ij} \geq 0 \qquad \forall i,j=1..n; i \neq j$

Comentario: N/A
"""

# ╔═╡ a0d7a86b-ad8d-4ad8-8d73-bc8f17940138
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
md"""
### Generación de Instancia
"""

# ╔═╡ 9ef5646a-9b56-4f72-a68d-43f4fc60d74a
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
begin
    c = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    n < 11 ? c : nothing
end

# ╔═╡ 249acea1-9778-4e03-9a7e-5636054c3f06
md"""
### Modelo en JuMP
"""

# ╔═╡ 3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
begin
    m = JuMP.Model()

    @variable(m, x[1:n, 1:n], Bin)
    @variable(m, f[1:n, 1:n] >= 0) # Gavish and Graves formulation

    @objective(m, Min, sum(c .* x))

    @constraint(m, r0[i in 1:n], x[i, i] == 0) # FIX
    @constraint(m, r1[i in 1:n], sum(x[i, :]) == 1)
    @constraint(m, r2[j in 1:n], sum(x[:, j]) == 1)

    # SEC Gavish and Graves formulation
    @constraint(m, r3, sum(f[1, :]) == n)
    @constraint(m, r4, sum(f[:, 1]) == 1)
    @constraint(m, r5[i in 1:n, j in 1:n], x[i, j] <= f[i, j])
    @constraint(m, r6[i in 1:n, j in 1:n], f[i, j] <= n * x[i, j])
    @constraint(m, r7[k in 2:n], sum(f[:, k]) == 1 + sum(f[k, :]))

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

# ╔═╡ 9b99a75d-00e0-4606-a160-3ceb14d15379
JuMP.value.(f)

# ╔═╡ db1b2c35-1196-485d-9d16-790cbfd9e1e7
fval = round.(Int, JuMP.value.(f))

# ╔═╡ ddd6790a-7597-4983-909d-41c11b24110c
md"""
### Solución gráfica
"""

# ╔═╡ c806912b-8bec-4d41-a7f3-0c7654a4ed53
begin
    tour = Int[1]
    while true
        push!(tour, argmax(xval[tour[end], :]))
        if tour[end] == 1
            break
        end
    end
    tour
end

# ╔═╡ 33c830a1-9fdc-430a-85c0-988723569738
sum(c[tour[i], tour[i+1]] for i in 1:n)

# ╔═╡ c2c05477-7c33-480f-b88f-c173487a4a4c
begin
    plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in 1:n
        annotate!(X[i], Y[i], text("$i", :top))
    end

    plot!(X[tour], Y[tour], color = :black)
end

# ╔═╡ Cell order:
# ╠═e7262522-ac65-11ec-0633-1d82420161db
# ╠═c22a2d8c-885b-4b52-a411-4eaf6d29c059
# ╠═791783f2-1f1f-44cc-bd7e-8425e0d18c2f
# ╠═a0d7a86b-ad8d-4ad8-8d73-bc8f17940138
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═9ef5646a-9b56-4f72-a68d-43f4fc60d74a
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═2edd33de-ebd4-484c-a45b-3fe7c6beb597
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╠═21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╠═9b99a75d-00e0-4606-a160-3ceb14d15379
# ╠═db1b2c35-1196-485d-9d16-790cbfd9e1e7
# ╠═ddd6790a-7597-4983-909d-41c11b24110c
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═33c830a1-9fdc-430a-85c0-988723569738
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═8b45efd6-0fec-4a1f-a5a6-2fa8d95e9a99

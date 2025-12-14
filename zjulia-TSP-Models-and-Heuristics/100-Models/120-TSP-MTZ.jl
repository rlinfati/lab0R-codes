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

# ╔═╡ 8bbdfc85-f462-4df1-8ebe-264fa14c8b6b
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ bb6b3c2d-a373-455b-941c-d8412440b612
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Miller-Tucker-Zemlin

Lectura recomendada:

```
Miller, C. E., Tucker, A. W., & Zemlin, R. A. (1960). Integer programming formulation of traveling salesman problems. Journal of the ACM (JACM), 7(4), 326-329.
```
"""

# ╔═╡ da486799-48c9-4045-8fda-0d086fd68725
md"""
## Modelo Matemático

Parametro $n$: Cantidad de ciudades.

Parametro $c_{ij}$: Costo (distancia, tiempo, etc.) de recorrer desde la ciudad $i$ a la ciudad $j$. Generalmente distancia euclidiana.

Variable $x_{ij}$ : Variable binaria que toma el valor 1 si el vendedor recorrer desde la ciudad $i$ a la ciudad $j$. 0 en otro caso (eoc).

Variable $u_{i}$ : Variable continua, etiqueta las ciudades con el orden de visita.

$Min \: z = \sum_{i=1}^{n} \sum_{\displaylines{j=1 \\ i \neq j}}^{n} c_{ij}x_{ij}$

$\sum_{\displaylines{i=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall j = 1..n$

$\sum_{\displaylines{j=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall i = 1..n$

$u_1 = 0$

$u_i \leq n-1 \qquad \forall i=2,..,n$

$u_i + 1 \leq u_j + n(1-x_{ij}) \qquad \forall i,j=1..n, i \neq j$ 

$x_{ij} \in \{0,1\} \qquad \forall i,j=1..n; i \neq j$

$u_i \geq 0 \qquad \forall i=1..n$

Comentario: La cota inferior (relajación lineal) es muy mala.
"""

# ╔═╡ ad1d41ae-e87e-4ffb-97d6-4b5a06ad4672
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
md"""
### Generación de Instancia
"""

# ╔═╡ d9ff363a-edc5-4525-9357-42246e030b83
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
    @variable(m, u[1:n] >= 0) # Miller-Tucker-Zemlin formulation

    @objective(m, Min, sum(c .* x))

    @constraint(m, r0[i in 1:n], x[i, i] == 0) # FIX
    @constraint(m, r1[i in 1:n], sum(x[i, :]) == 1)
    @constraint(m, r2[j in 1:n], sum(x[:, j]) == 1)

    # SEC Miller-Tucker-Zemlin formulation
    @constraint(m, r3, u[1] == 0)
    @constraint(m, r4[i in 2:n], u[i] <= n - 1)
    @constraint(m, r5[i in 1:n, j in 2:n], u[i] + 1 <= u[j] + n * (1 - x[i, j]))

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

# ╔═╡ 2f0a61c3-ac96-4571-985a-9bfd54b954bc
uval = JuMP.value.(u)

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
# ╠═bb6b3c2d-a373-455b-941c-d8412440b612
# ╠═da486799-48c9-4045-8fda-0d086fd68725
# ╠═ad1d41ae-e87e-4ffb-97d6-4b5a06ad4672
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═d9ff363a-edc5-4525-9357-42246e030b83
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═2edd33de-ebd4-484c-a45b-3fe7c6beb597
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╟─21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╠═2f0a61c3-ac96-4571-985a-9bfd54b954bc
# ╠═ddd6790a-7597-4983-909d-41c11b24110c
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═33c830a1-9fdc-430a-85c0-988723569738
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═8bbdfc85-f462-4df1-8ebe-264fa14c8b6b

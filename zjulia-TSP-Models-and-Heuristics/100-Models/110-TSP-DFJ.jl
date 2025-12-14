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
            Pkg.PackageSpec("Plots")
            Pkg.PackageSpec("Combinatorics")
            Pkg.PackageSpec("PlutoUI")
        ],
    )
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

# ╔═╡ e71a75f8-6ba3-4571-8b33-d68cc28071dd
using Combinatorics

# ╔═╡ 7e0dcaba-b114-4a2b-bcb7-fcc8b48ed363
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ 66bea84e-7574-494c-abcf-eb8d092e778c
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Dantzig-Fulkerson-Johnson

Lectura recomendada:

```
Dantzig, G., Fulkerson, R., & Johnson, S. (1954). Solution of a large-scale traveling-salesman problem. Journal of the operations research society of America, 2(4), 393-410.
```
"""

# ╔═╡ ed843770-c34d-4e98-b61c-b9e58da449f6
md"""
## Modelo Matemático

Parametro $n$: Cantidad de ciudades.

Parametro $c_{ij}$: Costo (distancia, tiempo, etc.) de recorrer desde la ciudad $i$ a la ciudad $j$. Generalmente distancia euclidiana.

Variable $x_{ij}$ : Variable binaria que toma el valor 1 si el vendedor recorrer desde la ciudad $i$ a la ciudad $j$. 0 en otro caso (eoc).

$Min \: z = \sum_{i=1}^{n} \sum_{\displaylines{j=1 \\ i \neq j}}^{n} c_{ij}x_{ij}$

$\sum_{\displaylines{i=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall j = 1..n$

$\sum_{\displaylines{j=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall i = 1..n$

$\sum_{i \in Q}^{n} \sum_{\displaylines{j \in Q \\ i \neq j}}^{n} x_{ij} \leq |Q|-1 \qquad \forall Q \subsetneqq \{1,..,n\}, |Q| \geq 2$

$x_{ij} \in \{0,1\} \qquad \forall i,j=1..n; i \neq j$

Comentario: Q es un subconjunto de ciudades, debe contener mínimo dos ciudades, pero no puede contenerlas todas. Existen MUCHOS Q, una cantidad exponencial.
"""

# ╔═╡ 9281043c-05d4-4709-8e18-20b609d39d38
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
md"""
### Generación de Instancia
"""

# ╔═╡ 8aa4ee2a-ef01-4f2f-96c4-1d393ad9c2b8
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

    @objective(m, Min, sum(c .* x))

    @constraint(m, r0[i in 1:n], x[i, i] == 0) # FIX
    @constraint(m, r1[i in 1:n], sum(x[i, :]) == 1)
    @constraint(m, r2[j in 1:n], sum(x[:, j]) == 1)

    # SEC Dantzig-Fulkerson-Johnson formulation
    for s in powerset(1:n, 2, n - 1)
        @constraint(m, sum(x[i, j] for i in s, j in s) <= length(s) - 1)
    end

    n < 11 ? JuMP.latex_formulation(m) : nothing
end

# ╔═╡ e5acaf7f-4789-4c36-a959-b21ac389f913
n < 11 ? collect(powerset(1:n, 2, n - 1)) : nothing

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

# ╔═╡ cf5bcc3c-3e2b-43b8-b85e-ad1cd9928ad0
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
# ╠═66bea84e-7574-494c-abcf-eb8d092e778c
# ╠═ed843770-c34d-4e98-b61c-b9e58da449f6
# ╠═9281043c-05d4-4709-8e18-20b609d39d38
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═e71a75f8-6ba3-4571-8b33-d68cc28071dd
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═8aa4ee2a-ef01-4f2f-96c4-1d393ad9c2b8
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═e5acaf7f-4789-4c36-a959-b21ac389f913
# ╠═2edd33de-ebd4-484c-a45b-3fe7c6beb597
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╟─21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╟─ddd6790a-7597-4983-909d-41c11b24110c
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═cf5bcc3c-3e2b-43b8-b85e-ad1cd9928ad0
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═7e0dcaba-b114-4a2b-bcb7-fcc8b48ed363

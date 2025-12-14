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

# ╔═╡ 12120148-bf09-4215-a58a-f30549cab3bd
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ 66bea84e-7574-494c-abcf-eb8d092e778c
md"""
# Traveling Salesman Problem (TSP)
Tambien conocido como:
- Problema del vendedor viajero
- Problema del vendedor ambulante
- Problema del agente de comercio
- Problema del viajante

## Definición
Dada una lista de ciudades y las distancias entre cada par de ellas, ¿cuál es la ruta más corta posible que visita cada ciudad exactamente una vez y al finalizar regresa a la ciudad origen?

[Más detalles en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

[The Traveling Salesman Problem: Postcards from the Edge of Impossibility](https://www.youtube.com/watch?v=5VjphFYQKj8)
"""

# ╔═╡ b0c8f001-55b0-4493-ba79-d8e569191938
md"""
## Modelo Matemático

Parametro $n$: Cantidad de ciudades.

Parametro $c_{ij}$: Costo (distancia, tiempo, etc.) de recorrer desde la ciudad $i$ a la ciudad $j$. Generalmente distancia euclidiana.

Variable $x_{ij}$ : Variable binaria que toma el valor 1 si el vendedor recorrer desde la ciudad $i$ a la ciudad $j$. 0 en otro caso (eoc).

$$\require{ams}$$

$Min \: z = \sum_{i=1}^{n} \sum_{\displaylines{j=1 \\ i \neq j}}^{n} c_{ij}x_{ij}$

$\sum_{\displaylines{i=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall j = 1..n$

$\sum_{\displaylines{j=1 \\ i \neq j}}^{n} x_{ij} = 1 \qquad \forall i = 1..n$

$Subtour \: Elimination \: Constraints.$

$x_{ij} \in \{0,1\} \qquad \forall i,j=1..n; i \neq j$

Comentario: Faltan las SEC
"""

# ╔═╡ 2c24cad5-188f-4ba7-85b0-06bb3a3135ef
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
md"""
### Generación de Instancia
"""

# ╔═╡ 5cd05e18-4162-4dd9-a6d6-d8ae6ed4f94a
n = 10

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

    n < 11 ? JuMP.latex_formulation(m) : nothing
end

# ╔═╡ 5f646768-ff04-426e-822e-25d6c0c9d58d
md"""
### Parametros del Solver y Solución
"""

# ╔═╡ 3c1bc4e6-05a2-490a-8c0e-321df5053a05
begin
    JuMP.set_optimizer(m, GLPK.Optimizer)
    JuMP.set_optimizer_attribute(m, "msg_lev", GLP_MSG_ALL)
    JuMP.set_optimizer_attribute(m, "tm_lim", 60 * 1000)
    JuMP.optimize!(m)
end

# ╔═╡ 532d2020-7dff-4f70-b6fc-9523040f843f
begin
    xval0 = JuMP.value.(x) .≈ 1.0
    xpar0 = [[i, j] for i in 1:n, j in 1:n if xval0[i, j]]
end

# ╔═╡ 88da88a7-7bae-4789-a2d7-6b491f41dd89
begin
    p0 = plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in 1:n
        annotate!(X[i], Y[i], text("$i", :top))
    end

    for i in xpar0
        plot!(X[i], Y[i], color = :black)
    end

    p0
end

# ╔═╡ 4a5b5694-74b4-4414-90c0-b8c727015590
md"""
### Eliminacion de subtours
"""

# ╔═╡ 2fd4955b-f0ae-443d-9ed9-f8c41ff8c391
subtour = [[1, 5], [1, 5, 10, 4, 8], [1, 5, 6, 2, 3, 9, 7, 8], [1, 5, 10, 4]]

# ╔═╡ b29e3b73-9c66-450f-bcd2-d8ab378a5576
begin
    for st in subtour
        @constraint(m, sum(x[i, j] for i in st, j in st) <= length(st) - 1)
    end
    JuMP.optimize!(m)
    xval = JuMP.value.(x) .≈ 1.0
    xpar = [[i, j] for i in 1:n, j in 1:n if xval[i, j]]
end

# ╔═╡ 1291457d-fdc9-4556-8927-a00c66582ce5
JuMP.solution_summary(m)

# ╔═╡ 3c678fae-3763-4235-afb8-1f6fc7b3c302
begin
    p = plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in 1:n
        annotate!(X[i], Y[i], text("$i", :top))
    end

    for i in xpar
        plot!(X[i], Y[i], color = :black)
    end

    p
end

# ╔═╡ 21cc515d-d381-445f-a23a-4bc75c81e38c
md"""
### Solución en forma de tour
"""

# ╔═╡ 43479c42-4e4e-4cb2-a462-cad721f1b18f
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

# ╔═╡ ea89fc90-dacd-4f90-852b-a9eb3bce0450
sum(c[tour[i], tour[i+1]] for i in 1:n)

# ╔═╡ 370e3cc2-a763-498f-88e2-a871ddd30158
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
# ╠═b0c8f001-55b0-4493-ba79-d8e569191938
# ╠═2c24cad5-188f-4ba7-85b0-06bb3a3135ef
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═5cd05e18-4162-4dd9-a6d6-d8ae6ed4f94a
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═5f646768-ff04-426e-822e-25d6c0c9d58d
# ╠═3c1bc4e6-05a2-490a-8c0e-321df5053a05
# ╠═532d2020-7dff-4f70-b6fc-9523040f843f
# ╠═88da88a7-7bae-4789-a2d7-6b491f41dd89
# ╠═4a5b5694-74b4-4414-90c0-b8c727015590
# ╠═2fd4955b-f0ae-443d-9ed9-f8c41ff8c391
# ╠═b29e3b73-9c66-450f-bcd2-d8ab378a5576
# ╠═1291457d-fdc9-4556-8927-a00c66582ce5
# ╠═3c678fae-3763-4235-afb8-1f6fc7b3c302
# ╠═21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═43479c42-4e4e-4cb2-a462-cad721f1b18f
# ╠═ea89fc90-dacd-4f90-852b-a9eb3bce0450
# ╠═370e3cc2-a763-498f-88e2-a871ddd30158
# ╠═12120148-bf09-4215-a58a-f30549cab3bd

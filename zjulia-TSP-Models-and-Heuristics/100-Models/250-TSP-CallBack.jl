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

# ╔═╡ 9a40869c-f9a7-4054-aee3-a3ca077df8b0
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ 66bea84e-7574-494c-abcf-eb8d092e778c
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Formulación Miller-Tucker-Zemlin

Lectura recomendada:

```
Miller, C. E., Tucker, A. W., & Zemlin, R. A. (1960). Integer programming formulation of traveling salesman problems. Journal of the ACM (JACM), 7(4), 326-329.
```
"""

# ╔═╡ 88a0fdad-83fd-4079-a193-aeaf7058c8b5
md"""
## Implementación con JuMP y GLPK
"""

# ╔═╡ 4e21d0be-afd3-4e02-84a4-f4a39a2d712f
md"""
### Generación de Instancia
"""

# ╔═╡ 3c822ca2-3e94-4e6b-bbc1-2a9bd599d05b
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

    # SEC-2 Dantzig-Fulkerson-Johnson formulation
    for a in 1:n, b in a+1:n
        s = [a, b]
        @constraint(m, sum(x[i, j] for i in s, j in s) <= length(s) - 1)
    end

    n < 11 ? JuMP.latex_formulation(m) : nothing
end

# ╔═╡ 7cfdcf29-7d9c-4122-8ec1-64dfd1c23e1a
"""
Función tspDist: dada una matriz de distancia y un tour, calcula la distancia total
"""
function tspDist(dist::Array{Float64,2}, tour::Array{Int,1})
    n, _ = size(dist)

    f = tour[1] == tour[end] && sort(tour[1:(end-1)]) == 1:n
    if f == false
        return +Inf
    end

    totalDist = sum(dist[tour[i], tour[i+1]] for i in 1:n)
    return totalDist
end

# ╔═╡ fff51859-5a07-4c3b-881c-06fe01211aeb
"""
Función myCallbackUserCut: genera SEC-3 Dantzig-Fulkerson-Johnson fraccionarios
"""
function myCallbackUserCut(cb_data)
    ret = callback_node_status(cb_data, m)
    @show "C", ret

    xval = JuMP.callback_value.(cb_data, m[:x])

    for a in 1:n, b in a+1:n, c in b+1:n
        s = [a, b, c]
        sums = sum(xval[i, j] for i in s, j in s)

        if sums <= length(s) - 1 + eps()
            continue
        end

        @show "c", sums, s
        con = @build_constraint(sum(m[:x][i, j] for i in s, j in s) <= length(s) - 1)
        ret = MOI.submit(m, MOI.UserCut(cb_data), con)
    end

    return
end

# ╔═╡ 702938a7-922b-46e5-b0ba-bf977a02ac83
"""
Función tspCWs: dada una matriz de distancia, calcula un tour factible.
"""
function tspCWs(dist::Array{Float64,2}, lamda::Float64 = 1.0, s::Int = 0)
    n, _ = size(dist)

    # Heuristica de los ahorros
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

# ╔═╡ 523d44ac-351f-4386-915a-b8eaa39f10e2
"""
Función tsp2Opt: dado un tour factible, y una una matriz de distancia, calcula un tour 2-OPT mejorado
"""
function tsp2Opt(tour::Array{Int,1}, dist::Array{Float64,2}, sf::Bool = false)
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

# ╔═╡ e858baeb-4a31-4101-a45e-ee3834eb4cf6
"""
Función mhLocalSearch1: dado un tour factible, y una una matriz de distancia, calcula un tour optimo local
"""
function mhLocalSearch1(tour::Array{Int,1}, dist::Array{Float64,2}, sf::Bool = false)
    ctour = copy(tour)
    busqueda0 = [tspDist(dist, ctour)] # actual ztour

    while true
        newtour = tsp2Opt(ctour, dist, sf)
        if newtour === nothing
            break
        end

        newtourz = tspDist(dist, newtour)

        if newtourz < busqueda0[end]
            push!(busqueda0, newtourz)
            ctour = newtour
        else
            break
        end
    end

    return ctour
end

# ╔═╡ 2edd33de-ebd4-484c-a45b-3fe7c6beb597
md"""
### Referencias "globales" para callbacks
"""

# ╔═╡ 388b624b-b9c1-4477-8726-a2e4c57a5566
begin
    global cbh_firstrun::Ref{Bool} = true
    global tourLz::Ref{Float64} = sum(maximum(c, dims = 2))
    global tourL::Ref{Array{Int}} = []
    nothing
end

# ╔═╡ 7b9ba7f4-2281-4d81-86bb-a5b71e10e887
"""
Función myCallbackLazyConstraint: genera SEC Dantzig-Fulkerson-Johnson enteras
"""
function myCallbackLazyConstraint(cb_data)
    ret = callback_node_status(cb_data, m)
    @show "L", ret
    if ret == JuMP.MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        return
    end

    xval = JuMP.callback_value.(cb_data, m[:x])

    aVisitar = collect(1:n)
    while !isempty(aVisitar)
        s = Int[aVisitar[1]]
        while true
            proxID = argmax(xval[s[end], aVisitar])

            if s[1] != aVisitar[proxID]
                push!(s, aVisitar[proxID])
                deleteat!(aVisitar, proxID)
            else
                deleteat!(aVisitar, proxID)
                break
            end
        end

        if length(s) == n
            ss = [s..., s[1]]
            zTSP = tspDist(c, ss)
            if zTSP < tourLz[]
                tourLz[] = zTSP
                tourL[] = ss
            end
            return
        end

        @show "l", length(s), s
        con = @build_constraint(sum(m[:x][i, j] for i in s, j in s) <= length(s) - 1)
        MOI.submit(m, MOI.LazyConstraint(cb_data), con)
    end
    return
end

# ╔═╡ f2f892e0-52ef-4945-a8b7-36e4fb15b44f
"""
Función myCallbackHeuristic: genera una solucion heuristica para el branch-and-bound
"""
function myCallbackHeuristic(cb_data)
    ret = callback_node_status(cb_data, m)
    @show "H", ret

    if cbh_firstrun[]
        tour = tspCWs(c)
        tour = mhLocalSearch1(tour, c)

        xval = zeros(n, n)
        for i in 1:n
            xval[tour[i], tour[i+1]] = 1.0
        end

        ret = MOI.submit(m, MOI.HeuristicSolution(cb_data), x[:], xval[:])
        println("** myCallbackHeuristic0 status = $(ret)")

        cbh_firstrun[] = false
    end

    if length(tourL[]) == n + 1
        tour = mhLocalSearch1(tourL[], c)
        tourz = tspDist(c, tour)

        if tourz < tourLz[]
            xval = zeros(n, n)
            for i in 1:n
                xval[tour[i], tour[i+1]] = 1.0
            end

            ret = MOI.submit(m, MOI.HeuristicSolution(cb_data), x[:], xval[:])
            println("** myCallbackHeuristicL status = $(ret)")
        end
        tourLz[] = tourz - eps(Float16)
        tourL[] = Int[]
    end

    return
end

# ╔═╡ dfcbe57c-9fd1-4f41-886d-d57581d4403b
md"""
### Parametros del Solver y Solución
"""

# ╔═╡ d1b599d9-da3a-4070-90dc-e63532951fd6
begin
    cbh_firstrun[] = true
    tourLz[] = sum(sum(maximum(c, dims = 2)))
    tourL[] = Int[]

    JuMP.set_optimizer(m, GLPK.Optimizer)
    JuMP.set_optimizer_attribute(m, "msg_lev", GLP_MSG_ALL)
    JuMP.set_optimizer_attribute(m, "tm_lim", 60 * 1000)

    JuMP.MOI.set(m, JuMP.MOI.LazyConstraintCallback(), myCallbackLazyConstraint)
    JuMP.MOI.set(m, JuMP.MOI.UserCutCallback(), myCallbackUserCut)
    JuMP.MOI.set(m, JuMP.MOI.HeuristicCallback(), myCallbackHeuristic)

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

# ╔═╡ dc506999-2c80-42ed-a2c2-84b6aad636d2
sum(c[tour[i], tour[i+1]] for i in 1:n)

# ╔═╡ c2c05477-7c33-480f-b88f-c173487a4a4c
begin
    p = plot(legend = false)
    scatter!(X, Y, color = :blue)

    for i in 1:n
        annotate!(X[i], Y[i], text("$i", :top))
    end

    plot!(X[tour], Y[tour], color = :black)
end

# ╔═╡ Cell order:
# ╠═e7262522-ac65-11ec-0633-1d82420161db
# ╠═66bea84e-7574-494c-abcf-eb8d092e778c
# ╠═88a0fdad-83fd-4079-a193-aeaf7058c8b5
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a28468d5-b6fd-42d1-856f-4b9a8198a4fa
# ╠═e6902e75-32f3-4c61-b564-24ce502f7025
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═4e21d0be-afd3-4e02-84a4-f4a39a2d712f
# ╠═3c822ca2-3e94-4e6b-bbc1-2a9bd599d05b
# ╠═f9cb6029-062b-436e-968b-9e8e4dd7b8f6
# ╠═72b74a0b-e001-445f-9c65-995c69854951
# ╠═c2a2b9e2-170f-44eb-be42-ea3219836e91
# ╠═249acea1-9778-4e03-9a7e-5636054c3f06
# ╠═3f77731c-9b9e-4ac2-bd91-cc7f8e528faf
# ╠═7cfdcf29-7d9c-4122-8ec1-64dfd1c23e1a
# ╠═7b9ba7f4-2281-4d81-86bb-a5b71e10e887
# ╠═fff51859-5a07-4c3b-881c-06fe01211aeb
# ╠═702938a7-922b-46e5-b0ba-bf977a02ac83
# ╠═523d44ac-351f-4386-915a-b8eaa39f10e2
# ╠═e858baeb-4a31-4101-a45e-ee3834eb4cf6
# ╠═f2f892e0-52ef-4945-a8b7-36e4fb15b44f
# ╠═2edd33de-ebd4-484c-a45b-3fe7c6beb597
# ╠═388b624b-b9c1-4477-8726-a2e4c57a5566
# ╠═dfcbe57c-9fd1-4f41-886d-d57581d4403b
# ╠═d1b599d9-da3a-4070-90dc-e63532951fd6
# ╠═fc925cdf-ffd8-45ad-a7a9-4c11228fac02
# ╟─21cc515d-d381-445f-a23a-4bc75c81e38c
# ╠═2fd0b79d-4740-45aa-a63d-95d2fb41ad81
# ╠═ddd6790a-7597-4983-909d-41c11b24110c
# ╠═c806912b-8bec-4d41-a7f3-0c7654a4ed53
# ╠═dc506999-2c80-42ed-a2c2-84b6aad636d2
# ╠═c2c05477-7c33-480f-b88f-c173487a4a4c
# ╠═9a40869c-f9a7-4054-aee3-a3ca077df8b0

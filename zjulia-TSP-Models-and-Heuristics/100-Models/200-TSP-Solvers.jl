### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ f6ca6701-c1f3-436a-9f7a-a1419d7a75b2
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    Pkg.add(
        [
            Pkg.PackageSpec("JuMP")
            Pkg.PackageSpec("GLPK")   # GPL
            Pkg.PackageSpec("HiGHS")  # MIT
            Pkg.PackageSpec("SCIP")   # ZIB
            Pkg.PackageSpec("COPT")   # Commercial
            Pkg.PackageSpec("Gurobi") # commercial
            #Pkg.PackageSpec("CPLEX")  # commercial
            Pkg.PackageSpec("PlutoUI")
        ],
    )
    Pkg.status()
end

# ╔═╡ e726254a-ac65-11ec-2f4c-1bead54e006a
using JuMP

# ╔═╡ a6bc0450-5aee-4216-a54c-0a64fcf114f1
using GLPK, HiGHS, SCIP

# ╔═╡ 22c0bbcf-0d00-4f0b-97ba-3a49df97cabb
using COPT, Gurobi #CPLEX

# ╔═╡ bb4218e4-06f0-482b-9ab6-55fedd51d429
using Random

# ╔═╡ 1f5153e3-3138-4864-86dc-1157212e2dc7
using Plots

# ╔═╡ 79445d8a-e0bb-496a-8e51-a59306c89e6d
begin
    using PlutoUI
    PlutoUI.TableOfContents()
end

# ╔═╡ e7262522-ac65-11ec-0633-1d82420161db
md"""
# Traveling Salesman Problem (TSP)
[TSP en Wikipedia](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## Open Source Solver

- [GLPK](http://www.gnu.org/software/glpk/) + [GLPK.jl](https://github.com/jump-dev/GLPK.jl)
- [HiGHS](https://github.com/ERGO-Code/HiGHS) + [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)
- [SCIP](https://scipopt.org/) + [SCIP.jl](https://github.com/scipopt/SCIP.jl)

## Closed Source Solver
- [COPT](https://www.shanshu.ai/copt) + [COPT.jl](https://github.com/COPT-Public/COPT.jl)
- [Gurobi](https://gurobi.com) + [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/) + [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)
"""

# ╔═╡ 22b69233-3007-4352-8551-0aa3af45f149
md"""
## Implementación con JuMP y varios solvers
"""

# ╔═╡ 65b95e72-2186-4a78-a0cf-49ebb4cfccf6
md"""
### Definicion de solvers y sus parametros
"""

# ╔═╡ 9d0c2f79-93f2-4bfe-a449-c8c2483a0016
function solverGLPK!(m::JuMP.Model)
    JuMP.set_optimizer(m, GLPK.Optimizer)
    JuMP.set_optimizer_attribute(m, "msg_lev", GLP_MSG_ALL)
    JuMP.set_optimizer_attribute(m, "tm_lim", 60 * 1000)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ 4457c6f1-bf75-4c47-814e-be4adda8a1a2
function solverHiGHS!(m::JuMP.Model)
    JuMP.set_optimizer(m, HiGHS.Optimizer)
    JuMP.set_optimizer_attribute(m, "log_to_console", true)
    JuMP.set_optimizer_attribute(m, "log_file", joinpath(@__DIR__, "logfile.txt"))
    JuMP.set_optimizer_attribute(m, "time_limit", 60.0)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ 5175db90-1bf4-44c1-9ff0-65f81a25b138
function solverSCIP!(m::JuMP.Model)
    JuMP.set_optimizer(m, SCIP.Optimizer)
    JuMP.set_optimizer_attribute(m, "display/verblevel", 3)
    JuMP.set_optimizer_attribute(m, "limits/time", 60)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ 5de47f3e-63f0-4a5f-98db-e7701e407f27
function solverCOPT!(m::JuMP.Model)
    JuMP.set_optimizer(m, COPT.Optimizer)
    JuMP.set_optimizer_attribute(m, "LogToConsole", 1)
    JuMP.set_optimizer_attribute(m, "TimeLimit", 60)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ 48e6f211-27a6-42a2-8c4a-4352863da1fa
function solverGurobi!(m::JuMP.Model)
    JuMP.set_optimizer(m, Gurobi.Optimizer)
    JuMP.set_optimizer_attribute(m, "LogToConsole", 1)
    JuMP.set_optimizer_attribute(m, "LogFile", joinpath(@__DIR__, "logfile.txt"))
    JuMP.set_optimizer_attribute(m, "TimeLimit", 60)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ 6db6ca58-bd42-40f9-9ff6-639991ed5615
function solverCPLEX!(m::JuMP.Model)
    JuMP.set_optimizer(m, CPLEX.Optimizer)
    JuMP.set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", 3)
    JuMP.set_optimizer_attribute(m, "CPX_PARAM_TILIM", 60)
    JuMP.optimize!(m)
    return nothing
end

# ╔═╡ cc07b58f-025b-4b30-8b65-cd8f7c8bea4b
md"""
### Generación de Instancia
"""

# ╔═╡ 7d522813-71ef-41b4-91ac-2e587a874ec7
function instance01(n::Int)
    Random.seed!(1234)
    X = rand(n) * 1_000.0
    Y = rand(n) * 1_000.0
    return X, Y
end

# ╔═╡ 88fbbca4-4695-4bf0-8633-869251b76e62
md"""
### Calculo de matriz de costos
"""

# ╔═╡ d4097f60-aad1-4a62-a90d-3c662eefd744
function processInsta(X::Vector{Float64}, Y::Vector{Float64})
    n = length(X)
    @assert n == length(Y)

    c = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return c
end

# ╔═╡ 83d8a649-13e3-47da-949c-08f99fc2113e
md"""
### Modelo en JuMP
"""

# ╔═╡ 7f63ea6c-02b9-4920-8b2b-997b3ca01cbb
function modelTSP(c::Matrix{Float64})
    n, n2 = size(c)
    @assert n == n2
    @assert sum(c[i, i] for i in 1:n) ≈ 0.0

    m = JuMP.Model()

    @variable(m, x[1:n, 1:n], Bin)
    @variable(m, u[1:n] >= 0)

    @objective(m, Min, sum(c .* x))

    @constraint(m, r0[i in 1:n], x[i, i] == 0) # FIX
    @constraint(m, r1[i in 1:n], sum(x[i, :]) == 1)
    @constraint(m, r2[j in 1:n], sum(x[:, j]) == 1)

    # SEC Miller-Tucker-Zemlin formulation
    @constraint(m, r3, u[1] == 0)
    @constraint(m, r4[i in 2:n], u[i] <= n - 1)
    @constraint(m, r5[i in 1:n, j in 2:n], u[i] + 1 <= u[j] + n * (1 - x[i, j]))

    return m
end

# ╔═╡ 5789a581-f3c5-464c-8c7e-623e57d5f5b1
md"""
### Extrae tour desde solucion de JuMP
"""

# ╔═╡ 83f4e638-1488-467c-ae8e-f62f583c7ab3
function solutionTSP(m::JuMP.Model)
    if JuMP.primal_status(m) != JuMP.MOI.FEASIBLE_POINT
        return Int[]
    end

    xval = round.(Bool, JuMP.value.(m[:x]))

    tour = Int[1]
    while true
        push!(tour, argmax(xval[tour[end], :]))
        if tour[end] == 1
            break
        end
    end
    return tour
end

# ╔═╡ 3ebd7c94-14a4-420c-9b1a-8e456853f924
md"""
### Dibuja solucion desde tour
"""

# ╔═╡ 594761f6-941d-45f0-8c75-c48acda46b03
function plotTSP(X::Vector{Float64}, Y::Vector{Float64}, t::Vector{Int})
    p = plot(legend = false)
    scatter!(X, Y, color = :blue)

    if length(t) == length(X) + 1
        plot!(X[t], Y[t], color = :red)
    end

    return p
end

# ╔═╡ ec9f21a7-0a6e-4149-8cd4-58bda41decbb
md"""
## Ejemplo
### Crea instancia y resuelve instancia
"""

# ╔═╡ 9e869925-3a23-420f-aba3-8a83995012d1
n = 7

# ╔═╡ 7389e858-88c7-45db-9123-875dea2d22e5
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    solverGLPK!(m)
    println(JuMP.solution_summary(m))
    t = solutionTSP(m)
    p = plotTSP(X, Y, t)
end

# ╔═╡ 955c2310-e458-4b81-ab79-bf7235c56762
md"""
### solverGLPK
"""

# ╔═╡ 681a44f4-4c59-4f4d-96b7-7ea6eb3b0279
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    solverGLPK!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ ccfa1a52-08b1-408a-9041-dfcc9ce8f562
md"""
### solverHiGHS
"""

# ╔═╡ 64e9f958-ef77-47b7-b1a8-134fe6d164bc
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    solverHiGHS!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ 7b233b26-baaa-4705-abe0-87d929b2eff0
md"""
### solverSCIP
"""

# ╔═╡ e73a2dce-458f-42aa-82fe-af92cb35e8ad
md"""
### `Error: no BLAS/LAPACK library loaded!`
```julia
import LinearAlgebra, OpenBLAS32_jll
LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
```
"""

# ╔═╡ 31772c46-dc04-479d-914b-949e62a42364
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    solverSCIP!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ 01ddc508-80f1-4bc6-9b9b-4462ad5caa77
md"""
### solverCOPT
"""

# ╔═╡ 9b7c53de-afe9-4dd0-8611-31fe960eb87c
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    solverCOPT!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ 4a99bd13-c585-430a-a0dc-e31446a7e9b9
md"""
### solverGurobi
"""

# ╔═╡ d961a213-cc5c-469f-9e45-05a835d7f154
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    #solverGurobi!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ 99cdbb7b-dbcc-4168-97c5-ef4d89093e7c
md"""
### solverCPLEX
"""

# ╔═╡ d87d9814-3fec-4bbf-bfff-c6037b2576f6
let
    X, Y = instance01(n)
    d = processInsta(X, Y)
    m = modelTSP(d)
    #solverCPLEX!(m)
    println(JuMP.solution_summary(m))
end

# ╔═╡ Cell order:
# ╠═f6ca6701-c1f3-436a-9f7a-a1419d7a75b2
# ╠═e7262522-ac65-11ec-0633-1d82420161db
# ╠═22b69233-3007-4352-8551-0aa3af45f149
# ╠═e726254a-ac65-11ec-2f4c-1bead54e006a
# ╠═a6bc0450-5aee-4216-a54c-0a64fcf114f1
# ╠═22c0bbcf-0d00-4f0b-97ba-3a49df97cabb
# ╠═bb4218e4-06f0-482b-9ab6-55fedd51d429
# ╠═1f5153e3-3138-4864-86dc-1157212e2dc7
# ╠═65b95e72-2186-4a78-a0cf-49ebb4cfccf6
# ╠═9d0c2f79-93f2-4bfe-a449-c8c2483a0016
# ╠═4457c6f1-bf75-4c47-814e-be4adda8a1a2
# ╠═5175db90-1bf4-44c1-9ff0-65f81a25b138
# ╠═5de47f3e-63f0-4a5f-98db-e7701e407f27
# ╠═48e6f211-27a6-42a2-8c4a-4352863da1fa
# ╠═6db6ca58-bd42-40f9-9ff6-639991ed5615
# ╠═cc07b58f-025b-4b30-8b65-cd8f7c8bea4b
# ╠═7d522813-71ef-41b4-91ac-2e587a874ec7
# ╠═88fbbca4-4695-4bf0-8633-869251b76e62
# ╠═d4097f60-aad1-4a62-a90d-3c662eefd744
# ╠═83d8a649-13e3-47da-949c-08f99fc2113e
# ╠═7f63ea6c-02b9-4920-8b2b-997b3ca01cbb
# ╠═5789a581-f3c5-464c-8c7e-623e57d5f5b1
# ╠═83f4e638-1488-467c-ae8e-f62f583c7ab3
# ╠═3ebd7c94-14a4-420c-9b1a-8e456853f924
# ╠═594761f6-941d-45f0-8c75-c48acda46b03
# ╠═ec9f21a7-0a6e-4149-8cd4-58bda41decbb
# ╠═9e869925-3a23-420f-aba3-8a83995012d1
# ╠═7389e858-88c7-45db-9123-875dea2d22e5
# ╠═955c2310-e458-4b81-ab79-bf7235c56762
# ╠═681a44f4-4c59-4f4d-96b7-7ea6eb3b0279
# ╠═ccfa1a52-08b1-408a-9041-dfcc9ce8f562
# ╠═64e9f958-ef77-47b7-b1a8-134fe6d164bc
# ╠═7b233b26-baaa-4705-abe0-87d929b2eff0
# ╠═e73a2dce-458f-42aa-82fe-af92cb35e8ad
# ╠═31772c46-dc04-479d-914b-949e62a42364
# ╠═01ddc508-80f1-4bc6-9b9b-4462ad5caa77
# ╠═9b7c53de-afe9-4dd0-8611-31fe960eb87c
# ╠═4a99bd13-c585-430a-a0dc-e31446a7e9b9
# ╠═d961a213-cc5c-469f-9e45-05a835d7f154
# ╠═99cdbb7b-dbcc-4168-97c5-ef4d89093e7c
# ╠═d87d9814-3fec-4bbf-bfff-c6037b2576f6
# ╠═79445d8a-e0bb-496a-8e51-a59306c89e6d

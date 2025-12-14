### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 70964b64-bba1-11ec-2d9e-c9b4a3f47e09
begin
    import Pkg
    Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
    Pkg.activate()
    Pkg.add([
        Pkg.PackageSpec("DataFrames")
        Pkg.PackageSpec("JSON")
        Pkg.PackageSpec("XLSX")
    ])
    Pkg.status()
end

# ╔═╡ d3ab311b-5035-4083-92b1-905cfd2b8280
using DataFrames

# ╔═╡ b5d830df-465d-4990-8dd4-4396e90d4106
using JSON

# ╔═╡ 9a0c4276-d2da-426a-8844-c14bb6df4dc1
using XLSX

# ╔═╡ 154fae4f-799c-4868-9d05-7a732d5e8e68
using Downloads

# ╔═╡ 730585c8-b346-4a98-8cfc-a1de7592e13e
using Dates

# ╔═╡ 70964b9e-bba1-11ec-3155-0b491028e886
XLSX.openxlsx("elmer-" * Dates.format(Dates.today(), "yyyy-mm-dd") * ".xlsx", mode = "w") do xf
    url1 = "http://www.elmercurio.com/inversiones/json/jsonTablaRankingFFMM.aspx"
    url2 = "http://www.elmercurio.com/inversiones/json/jsonTablaFull.aspx"
    r_ffmm = Downloads.download(url1)
    s_ffmm = read(r_ffmm, String)
    j_ffmm = JSON.parse(s_ffmm)
    df_ffmm = vcat(DataFrame.(j_ffmm["rows"])...)

    sxf = XLSX.addsheet!(xf, "ID-0")
    XLSX.writetable!(sxf, collect(DataFrames.eachcol(df_ffmm)), DataFrames.names(df_ffmm))

    df_ffmm_cols1 = [:Categ, :VarParIndUmes, :InvNetaUmes, :VarParIndU3meses, :InvNetaU3meses, :Patrimonio]
    df_ffmm1 = df_ffmm[!, df_ffmm_cols1]
    sxf = XLSX.addsheet!(xf, "ID-0 Par")
    XLSX.writetable!(sxf, collect(DataFrames.eachcol(df_ffmm1)), DataFrames.names(df_ffmm1))

    df_ffmm_cols2 = [:Categ, :RtbUltMes, :RtbUlt3Meses, :RtbUltAnio, :RtbUlt5anios]
    df_ffmm2 = df_ffmm[!, df_ffmm_cols2]
    sxf = XLSX.addsheet!(xf, "ID-0 Rt")
    XLSX.writetable!(sxf, collect(DataFrames.eachcol(df_ffmm2)), DataFrames.names(df_ffmm2))

    the_ffmm = DataFrame()
    for r in eachrow(df_ffmm)
        @show r[:Id], r[:Categ]
        r_ffmms = Downloads.download(url2 * "?idcategoria=" * r["Id"])
        s_ffmms = read(r_ffmms, String)
        j_ffmms = JSON.parse(s_ffmms)
        df = vcat(DataFrame.(j_ffmms["rows"])...)
        the_ffmm = vcat(the_ffmm, df)

        sxf = XLSX.addsheet!(xf, "ID-" * r[:Id])
        XLSX.writetable!(sxf, collect(DataFrames.eachcol(df)), DataFrames.names(df))
    end

    ffmmagps = unique(the_ffmm[!, :adm])
    for ffmmagp in ffmmagps
        @show ffmmagp
        df_rows = [v == ffmmagp for v in the_ffmm[!, :adm]]
        df2 = the_ffmm[df_rows, :]

        sxf = XLSX.addsheet!(xf, "AGP-" * ffmmagp)
        XLSX.writetable!(sxf, collect(DataFrames.eachcol(df2)), DataFrames.names(df2))
    end

    df_run = []

    # SANTANDER
    ## Future Wealth YOU
    ## Acciones Globales ESG YOU
    ## Go Acciones Selectas USA YOU
    push!(df_run, "10057-YOU",  "8090-YOU", "8488-YOU")
    ## Ahorro Mediano Plazo YOU
    ## Renta Corto Plazo YOU
    ## Money Market YOU
    push!(df_run, "9222-YOU", "8615-YOU", "8057-YOU")
    ## Gestión Activa Agresiva YOU
    ## Gestión Activa Equilibrio YOU
    ## Gestión Activa Prudente YOU
    push!(df_run, "9648-YOU", "9646-YOU", "9649-YOU")

    # BCI
    ## Cartera Dinámica Activa CLASI
    ## Cartera Dinámica Balanceada CLASI
    ## Cartera Dinámica Conservadora CLASI
    push!(df_run, "8640-CLASI", "8639-CLASI", "8638-CLASI")
    ## Cartera Dinámica Chile CLASI
    ## Cartera Dinámica Ahorro CLASI
    ## Cartera Dinámica Corto Plazo CLASI
    push!(df_run, "9511-CLASI", "9228-CLASI", "8731-CLASI")
    ## Cartera Patrimonial Activa INVER
    ## Cartera Patrimonial Balanceada INVER
    ## Cartera Patrimonial Conservadora INVER
    push!(df_run, "9060-INVER", "9062-INVER", "9063-INVER")
    ## Cartera Patrimonial Ahorro INVER
    ## Cartera Patrimonial Corto Plazo INVER
    push!(df_run, "9061-INVER", "8976-INVER")
    ## Mach Arriesgado DIGITAL
    ## Mach Moderado DIGITAL
    ## Mach Conservador DIGITAL
    push!(df_run, "10578-DIGITAL", "10579-DIGITAL", "10580-DIGITAL")
    ## Ahorro Digital DIGITAL
    ## Mach DIGITAL
    push!(df_run, "10478-DIGITAL", "10358-DIGITAL")
    ## Global Titan CLASI
    push!(df_run, "8710-CLASI")

    # BANCHILE
    ## Agresivo DIGITAL
    ## Moderado DIGITAL
    ## Conservador DIGITAL
    push!(df_run, "10058-DIGITAL", "10060-DIGITAL", "10059-DIGITAL")
    ## Renta Futura DIGITAL
    ## Renta Corto Plazo DIGITAL
    ## Disponible DIGITAL
    push!(df_run, "8357-DIGITAL", "8274-DIGITAL", "8052-DIGITAL")
    ## BlackRock ESG DIGITAL
    ## Global Accionario DIGITAL
    push!(df_run, "10173-DIGITAL", "8088-DIGITAL")
    ## Portafolio Agresivo Largo Plazo L
    ## Portafolio Moderado Largo Plazo L
    ## Portafolio Conservador Largo Plazo L
    push!(df_run, "8448-L", "9043-L", "8377-L")
    ## Portafolio Ahorro Corto Plazo L
    ## Portafolio Retorno Mediano Plazo L
    push!(df_run, "10519-L", "10520-L")
    ## Capital Empresarial A
    push!(df_run, "9022-A")

    # SCOTIA
    ## Portafolio Más Arriesgado WEB
    ## Portafolio Arriesgado WEB
    ## Portafolio Moderado WEB
    ## Portafolio Conservador WEB
    ## Portafolio Más Conservador WEB
    push!(df_run, "8740-WEB", "8116-WEB", "8741-WEB", "8886-WEB", "8304-WEB")
    ## Acciones Chile WEB
    ## Acciones Chile Mid Cap WEB
    push!(df_run, "8289-WEB", "8857-WEB")
    ## Acciones Sustentables Global WEB
    ## Acciones USA WEB
    ## Acciones Europa WEB
    ## Renta Variable Latam WEB
    push!(df_run, "8822-WEB", "8480-WEB", "8484-WEB", "8482-WEB")
    ## Deuda Largo Plazo UF WEB
    ## Deuda Mediano Plazo UF WEB
    ## Deuda Corto Plazo UF WEB
    push!(df_run, "9154-WEB", "8106-WEB", "8991-WEB")
    ## Deuda Mediano Plazo WEB
    ## Deuda Chile Flexible WEB 
    push!(df_run, "8292-WEB", "9021-WEB")
    ## Liquidez WEB
    push!(df_run, "8187-WEB")
    ## Real Estate Global WEB
    push!(df_run, "8604-WEB")

    # ITAU
    ## Mi Cartera Lanzada SIMPLE
    ## Mi Cartera Aventurera SIMPLE
    ## Mi Cartera Exploradora SIMPLE
    ## Mi Cartera Tranqui SIMPLE
    push!(df_run, "10064-SIMPLE", "10063-SIMPLE", "10021-SIMPLE", "10020-SIMPLE")
    ## Gestionado Agresivo F1
    ## Gestionado Moderado F1
    ## Gestionado Conservador F1
    push!(df_run, "8993-F1", "8992-F1", "8994-F1")
    ## ETF It Now S&P IPSA UNICA
    ## ETF It Now S&P IPSA ESG UNICA
    ## ETF It Now Dividendo UNICA
    push!(df_run, "9019-UNICA", "10068-UNICA", "9685-UNICA")

    # FINTUAL
    ## Risky Norris A
    ## Moderate Pitt A
    ## Conservative Clooney A
    ## Very Conservative Streep A
    push!(df_run, "9570-A", "9569-A", "9568-A", "9730-A")

    df_cols = [Symbol("Rentb1 mes"), :Rentb3m, :RentbY, :Rentb12m]
    push!(df_cols, :FondoFull, :adm, :Run, :TAC, :Fondo)
    push!(df_cols, :invNet1m, :invNet1y, :varpar1m, :varpar1Y, :patrim, :par)

    df_rows = [v in df_run for v in the_ffmm[!, :FondoFull]]
    df2 = the_ffmm[df_rows, df_cols]

    sxf = XLSX.addsheet!(xf, "FFMM - RL")
    return XLSX.writetable!(sxf, collect(DataFrames.eachcol(df2)), DataFrames.names(df2))
end

# ╔═╡ Cell order:
# ╠═70964b64-bba1-11ec-2d9e-c9b4a3f47e09
# ╠═d3ab311b-5035-4083-92b1-905cfd2b8280
# ╠═b5d830df-465d-4990-8dd4-4396e90d4106
# ╠═9a0c4276-d2da-426a-8844-c14bb6df4dc1
# ╠═154fae4f-799c-4868-9d05-7a732d5e8e68
# ╠═730585c8-b346-4a98-8cfc-a1de7592e13e
# ╠═70964b9e-bba1-11ec-3155-0b491028e886

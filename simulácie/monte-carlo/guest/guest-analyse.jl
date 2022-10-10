function analyse(path)
    file = open(path, "r")
    electrons = JSON.parse(file)
    close(file)

    Tls = [Tleave(electron[1][end], electron[3][end], 0.511, 1.602e-19, 4000 / (40 * 10e-6), 40 * 10e-6, electron[6][end-1]) for electron in electrons]
    TlsHist = fit(Histogram, Tls, 0:1:10e3)
    TlsHP = scatter(x=TlsHist.edges[1], y=TlsHist.weights)

    T₀s = vcat(getindex.(electrons, 1)...)
    T₀sHist = fit(Histogram, T₀s, 0:100:10e3)
    T₀sHP = scatter(x=T₀sHist.edges[1], y=T₀sHist.weights)

    Tcs = filter(x -> x != -1, vcat(getindex.(electrons, 2)...))
    TcsHist = fit(Histogram, Tcs, 0:100:10e3)
    TcsHP = scatter(x=TcsHist.edges[1], y=TcsHist.weights)

    ϕs = vcat(getindex.(electrons, 3)...)
    ϕsHist = fit(Histogram, ϕs, 0:1e-3:maximum(ϕs))
    ϕsHP = scatter(x=ϕsHist.edges[1], y=ϕsHist.weights)

    θs = vcat(getindex.(electrons, 4)...)
    θsHist = fit(Histogram, θs, 1:1e-3:1.4)
    θsHP = scatter(x=θsHist.edges[1], y=θsHist.weights)

    las = rad2deg.([leaveAng(electron[1][end], electron[3][end], 0.511, 1.602e-19, 4000 / (40 * 10e-6), 40 * 10e-6, electron[6][end-1]) for electron in electrons])
    lasHist = fit(Histogram, las, rad2deg(1):rad2deg(1e-3):maximum(las))
    lasHP = scatter(x=lasHist.edges[1], y=lasHist.weights)

    diffs = [diff(electron[end-1]) for electron in electrons]

    dbc = vcat(diffs...)
    dbcHist = fit(Histogram, dbc, 1e-5:1e-7:maximum(dbc))
    dbcHP = scatter(x=dbcHist.edges[1], y=dbcHist.weights)

    diffVsPos = [scatter(x = electron[6][2:end-1], y = diff) for (electron, diff) in zip(electrons, diffs)]

    electrons, TlsHP, T₀sHP, TcsHP, ϕsHP, θsHP, lasHP, dbcHP, diffVsPos
end
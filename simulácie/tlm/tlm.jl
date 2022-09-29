using Trapz
using CairoMakie
using DelimitedFiles
using Polynomials

CairoMakie.activate!(type = "svg")

#= layout of matrices: 
    ---> t
    |
    |
    v x
=#

# redefine log function to not throw error when argument is negative, but return NaN
function ntlog(arg)
    if arg < 0
        @warn "Negative logarithm argument: $arg"
        return NaN
    end

    log(arg)
end

# MCP-PMT parameters:

#d = 1e-5        # m; channel diameter
#A = 6.45e-6     # m²; single anode area

#S = π*(d/2)^2   # channel cross section area
#Nc = A/S        # approx. number of channels

Vs = 2400       # V
Vm = 1100       # V
C = 2.86e-13    # F
R = 6.09e8      # Ω
Lmcp = 1.5e-3   # m
G = log(Vs/Vm)
k = 1.16e4      # m⁻¹

Qs = Vs*C
L = k*Lmcp
vs = Vs/L 

# space and time discretization
tsize = 100
xsize = 100
sizes = (xsize, tsize)
tind = 1:tsize
xind = 1:xsize

xcoords = range(0, L, length = xsize)
times = range(0, 5e-8, length = tsize)

ψ₀ = fill(0, xsize)

# prepare the for loop
ϵs = [range(0.1, 1, length = 10); fill(1, 20)]
N = length(ϵs)
ϵind = 1:N

# input signal
iₗ = Polynomial(readdlm("output-signal-poly.csv")[:, 1]).(times)
#i₀ = hcat([[zeros(100); fill(1e-11, 100)] for i in 1:10]...)[1:tsize]

#= series of gauss shaped peaks
bell_f(a, b, c) = x -> a*exp(-((x-b)^2)/c)
bells = [bell_f(1e-11, b, 1e-10) for b in 0:1e-4:1e-3]
i₀f(t) = sum(map(bf -> bf(t), bells))
i₀ = i₀f.(times)=#

function calcGain(i₀)
    # input charge as function of time
    #Q₀ = [trapz(times[1:t], exp.(times[1:t] ./ (R*C)) .* i₀[1:t]) for t in tind]

    # initial gain
    g₀ = [exp(G*xcoords[x]) for x in xind, t in tind]

    gs = Array{Float64, 3}(undef, xsize, tsize, N+1)
    gs[:, :, 1] = g₀

    QQss = Array{Float64, 3}(undef, xsize, tsize, N)
    QwQss = Array{Float64, 2}(undef, tsize, N)
    ψs = Array{Float64, 3}(undef, xsize, tsize, N)

    # calculate the gain
    for i in ϵind
        g = gs[:, :, i]

        QQs = [trapz(times[1:t], exp.(times[1:t] ./(R*C)) .* i₀[1:t] .* g[x, 1:t])/Qs for x in xind, t in tind]
        QQss[:, :, i] = QQs

        QwQs = [trapz(xcoords, QQs[:, t])/L for t in tind]
        QwQss[:, i] = QwQs

        ψ = [(ψ₀[x] + (QwQs[t] - QQs[x, t])*ϵs[i])*exp(-times[t]/(R*C)) for x in xind, t in tind]
        ψs[:, :, i] = ψ

        new_g = [exp(G*xcoords[x] + trapz(xcoords[1:x], ntlog.(1 .+ ψ[1:x, t]))) for x in xind, t in tind]
        gs[:, :, i+1] = new_g
    end

    gs, QQss, QwQss, ψs
end

function backwardGain()
    h₀ = [exp(G*(xcoords[x] - L)) for x in xind, t in tind]

    hs = Array{Float64, 3}(undef, xsize, tsize, N+1)
    hs[:, :, 1] = h₀

    QQss = Array{Float64, 3}(undef, xsize, tsize, N)
    QwQss = Array{Float64, 2}(undef, tsize, N)
    ψs = Array{Float64, 3}(undef, xsize, tsize, N)

    for i in ϵind
        h = hs[:, :, i]

        QQs = [trapz(times[1:t], iₗ[1:t] .* h[x, 1:t])/Qs for x in xind, t in tind]
        QQss[:, :, i] = QQs

        QwQs = [trapz(xcoords, QQs[:, t])/L for t in tind]
        QwQss[:, i] = QwQs

        ψ = [ψ₀[x] + (QwQs[t] - QQs[x, t])*ϵs[i] for x in xind, t in tind]
        ψs[:, :, i] = ψ

        new_h = [exp(G*(xcoords[x]-L) - trapz(xcoords[x:end], ntlog.(1 .+ ψ[x:end, t]))) for x in xind, t in tind]
        hs[:, :, i+1] = new_h
    end

    hs, QQss, QwQss, ψs
end

hs = backwardGain()[1]
i₀ = iₗ .* hs[1, :, end]

gs, QQss, QwQss, ψs = calcGain(i₀)

# calculate the signal current i(x, t)
g_final = gs[:, :, end]
i = [i₀[t] * g_final[x, t] for x in xind, t in tind]

# calculate the excess volateg Φ(x, t)
ψ_final = ψs[:, :, end]
Φ = [vs*trapz(xcoords[1:x], ψ_final[1:x, t]) for x in xind, t in tind]

# calculate the channel voltage V(x, t)
V = [vs*xcoords[x] + Φ[x, t] for x in xind, t in tind]


function dumpSteps(cartInd, gs, QQss, QwQss, ψs)
    dir = "dump"
    mkpath(dir)

    writedlm("$dir/g", gs[cartInd...])
    writedlm("$dir/QQs", QQss[cartInd...])
    writedlm("$dir/QwQs", QwQss[cartInd[2], cartInd[3]])
    writedlm("$dir/psi", ψs[cartInd...])
end

function plotGains(n, xinds, tinds, legend=true)
    plt = Figure(resolution = (600, 460))
    ax = Axis(
        plt[1, 1],
        yscale = log10,
        xlabel = "z [mm]",
        ylabel = "g(z, t)",
        spinewidth = 3,
        xgridwidth = 3,
        ygridwidth = 3,
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 25,
        yticklabelsize = 25
        )

    for t in tinds
        lines!(ax, xcoords[xinds] ./k .* 1e3, gs[xinds, t, n], label = "t = $(round(times[t] * 1e9, sigdigits=2)) ns", linewidth=3)
    end

    legend && axislegend(ax, position=:rb, labelsize=20)
    plt
end

function plotPulse(xinds, tinds, z, legend = true)
    plt = Figure(resolution = (600, 460))
    ax = Axis(plt[1, 1],
        xlabel = "t [ns]",
        ylabel = "i($z, t) [A]",
        spinewidth = 3,
        xgridwidth = 3,
        ygridwidth = 3,
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 25,
        yticklabelsize = 25
        )

    for x in xinds
        lines!(ax, times[tinds] .* 1e9, i[x, tinds], label = "z = $(round(xcoords[x]/k *1e3, sigdigits=2)) mm", linewidth=3)
    end

    legend && axislegend(ax, position=:rb)
    plt
end
using CSV
using DataFrames
using CairoMakie
using Latexify

set_theme!(Theme(
    spinewidth = 3,
    xgridwidth = 3,
    ygridwidth = 3,
    xlabelsize = 25,
    ylabelsize = 25,
    xticklabelsize = 25,
    yticklabelsize = 25
))

c = 299792458      # m/s
h = 6.62607015e-34 # Js
λ = 405e-9         # m
τ = 2.2e-9         # s

N(P, f) = P*λ/h/c/f

data = DataFrame(CSV.File("PvsClockF.csv", delim = ';'))
data[!, :N] = N.(data.P, data.f)

exp_data = DataFrame("f [MHz]" => data.f .* 1e-6, "P [\$\\mu\$W]" => data.P .* 1e6, "N" => Int.(round.(data.N, sigdigits = 2)))

plt = Figure()
ax = Axis(plt[1, 1], xlabel = "Clock frequency [MHz]", ylabel = "Number of photons in one pulse", xscale = log10)
scatterlines!(ax, data.f .* 1e-6, data.N, linestyle = :dot, linewidth = 3)

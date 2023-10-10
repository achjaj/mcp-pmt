using CairoMakie
using DelimitedFiles
using JSON

set_theme!(Theme(
    spinewidth = 3,
    xgridwidth = 3,
    ygridwidth = 3,
    xlabelsize = 25,
    ylabelsize = 25,
    xticklabelsize = 25,
    yticklabelsize = 25
))

cd("output/") # change to directory with output data
xcoords = readdlm("xcoords.csv")[:, 1]
times = readdlm("times.csv")[:, 1]
params = JSON.parsefile("parameters.json", dicttype = Dict{Symbol, Float64})

function plotG(xind, tind, step, legend = true; xtoz = true)
    data = readdlm("gs/$step.csv")
    cnvf = xtoz ? params[:k] : 1

    plt = Figure()
    ax = Axis(
            plt[1, 1],
            yscale = log10,
            xlabel = (xtoz ? "z [mm]" : "x"),
            ylabel = "g($(xtoz ? "z" : "x"), t)"
        )

    for t in tind
        lines!(ax, xcoords[xind] ./ cnvf, data[xind, t], label = "t = $(round(times[t] * 1e9, sigdigits=2)) ns")
    end

    legend && axislegend(ax, position=:rb)

    plt
end
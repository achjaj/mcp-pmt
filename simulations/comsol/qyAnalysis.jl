using Serialization
using PlotlyJS

times = Float64[]
qys = Vector{Vector{Float64}}()
indicies = Vector{Vector{Int}}()

histograms = GenericTrace{Dict{Symbol, Any}}[]
titles = ("Position along the channel [m]", "Count")

function parseComsol(comsol_file::String)
    global times
    global qys
    global indicies

    open(comsol_file, "r") do io
        ln = 9
        # skip 9 lines
        for i in 1:9
            readline(io)
        end

        while !eof(io)
            line = readline(io)
            ln += 1

            line[1] != '%' && throw(ErrorException("Line $ln: Unexpected token: '$(line[1])'"))

            push!(times, parse(Float64, split(line, "t=")[end]))

            line = readline(io)
            ln += 1

            qy = parse.(Float64, split(line, "\t"))
            basei = collect(eachindex(qy))

            nani = findall(x -> isnan(x), qy)
            deleteat!(qy, nani)
            deleteat!(basei, nani)

            push!(qys, qy)
            push!(indicies, basei)
        end
    end
end

function save(path::String)
    global times
    global qys
    global histograms
    global indicies

    serialize(path, (times, qys, histograms, indicies))
end

function load(path::String)
    global times
    global qys
    global histograms
    global indicies

    times, qys, histograms, indicies = deserialize(path)

    !(length(times) == length(qys) == length(histograms) == length(indicies)) && throw(ErrorException("Saved data have different lengths!"))
end

function histos(xbins)
    global qys
    global histograms

    for qy in qys
        push!(histograms, histogram(x = qy, autobinx = false, xbins = xbins))
    end
end

function exportHistograms(outDir::String; ext::String = "png")
    global times
    global histograms
    global titles

    !isdir(outDir) && mkdir(outDir)
    !isempty(readdir(outDir)) && throw(ErrorException("'$outDir' is not empty!"))

    for (i, hist) in enumerate(histograms)
        label = "$(times[i])s"
        p = plot(hist, Layout(xaxis_title = titles[1], yaxis_title = titles[2], title = "t = $label"))
        
        savefig(p, joinpath(outDir, "t$label.$ext"))
    end
end

function showHist(i::Int)
    global times
    global histograms
    global titles

    plot(histograms[i], Layout(xaxis_title = titles[1], yaxis_title = titles[2], title = "t = $(times[i])s"))
end
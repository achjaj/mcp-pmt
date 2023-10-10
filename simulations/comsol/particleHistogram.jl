using Plots
using LaTeXStrings
include("parseExport.jl")

channelLen = 4.2e-4 # m
binEdges = range(0, stop = channelLen, length = 22)[2:end-1]
noBins = length(binEdges)
binWidth = round(binEdges[2] - binEdges[1], sigdigits = 1)

particles = parseExport("pivi/particle-pivi-1e4.txt")
timeIdxs = eachindex(particles[1].times)

Nhistograms = Vector{Vector{Int}}()
vHistograms = Vector{Vector{Vector{Float64}}}() # array of "histograms"; each bin of a "histogram" is an array which contains pidx of particles belonging to the bin

speedsInTime = Vector{Vector{Float64}}()
yInTime = Vector{Vector{Float64}}()

for timeIdx in timeIdxs
    Nhistogram = zeros(length(binEdges))
    vHistogram = [Float64[] for edge in binEdges]
    
    speeds = Float64[]
    ys = Float64[]

    for particle in particles
        qy = particle.features["qy"][timeIdx]
        V = particle.features["cpt.V"][timeIdx]

        if !isnan(qy)
            push!(speeds, V)
            push!(ys, qy)

            for (binIdx, edge) in enumerate(binEdges)
                if qy <= edge
                    Nhistogram[binIdx] += 1
                    push!(vHistogram[binIdx], particle.features["cpt.V"][timeIdx])
                    break
                end
            end
        end
    end

    push!(Nhistograms, Nhistogram)
    push!(vHistograms, vHistogram)

    push!(speedsInTime, speeds)
    push!(yInTime, ys)
end

maxN = maximum(maximum.(Nhistograms))
NhistAnim = @animate for (hi, hist) in enumerate(Nhistograms)
    plot(binEdges .* 1e3, hist, ylims = (0, maxN), xlabel = "Position along the channel [mm]", ylabel = "Number of particles", legend = false, linetype = :steppre, grid = false)
    annotate!(0.35, 2/3*maxN, text("t = ~$(round(particles[1].times[hi]*1e12, digits = 1))ps\nNumber of bins: $noBins\nBin width: ~$(binWidth*1e3)mm", pointsize = 12))
    #annotate!(0.35, 195, )
end
gif(NhistAnim, "/tmp/NhistAnim.gif", fps = 15)

avrgVHistograms = Vector{Vector{Float64}}()
for vHistogram in vHistograms
    avrgVHistogram = Float64[]

    for bin in vHistogram
        push!(avrgVHistogram, isempty(bin) ? 0 : sum(bin)/length(bin))
    end

    push!(avrgVHistograms, avrgVHistogram)
end

maxAvrgV = maximum(maximum.(avrgVHistograms)) *1e-7
avrgVHistAnim = @animate for (hi, hist) in enumerate(avrgVHistograms)
    plot(binEdges .* 1e3, hist .* 1e-7, xlabel = "Position along the channel [mm]", ylabel = "Average speed [m/s x 10⁷]", ylims = (0, maxAvrgV), legend = false, linetype = :steppre, grid = false)
    #annotate!(0.35, 1e7, "")
    annotate!(0.35, maxAvrgV/3, text("t = ~$(round(particles[1].times[hi]*1e12, digits = 1))ps\nNumber of bins: $noBins\nBin width: ~$(binWidth*1e3)mm", pointsize = 12))
end
gif(avrgVHistAnim, "/tmp/avrgVHistAnim.gif", fps = 15)

avrgVsqrHistograms = Vector{Vector{Float64}}()
for vHistogram in vHistograms
    avrgVsqrHistogram = Float64[]

    for bin in vHistogram
        push!(avrgVsqrHistogram, isempty(bin) ? 0 : sum(bin .^ 2)/length(bin))
    end

    push!(avrgVsqrHistograms, avrgVsqrHistogram)
end

maxAvrgVsqrt = maximum(maximum.(avrgVsqrHistograms))
avrgVsqrHistAnim = @animate for (hi, hist) in enumerate(avrgVsqrHistograms)
    plot(binEdges .* 1e3, hist, xlabel = "Position along the channel [mm]", ylabel = L"\langle v^2 \rangle \quad \textrm{[m^2/s^2]}", title = L"\textbf{Distribution\ of\ \langle v^2 \rangle\ along\ the\ channel}", ylims = (0, maxAvrgVsqrt), legend = false, linetype = :steppre, grid = false)
    #annotate!(0.35, 1e7^2, "")
    annotate!(0.35, maxAvrgVsqrt/3, text("t = ~$(round(particles[1].times[hi]*1e12, digits = 1))ps\nNumber of bins: $noBins\nBin width: ~$(binWidth*1e3)mm", pointsize = 12))
end
gif(avrgVsqrHistAnim, "/tmp/avrgVsqrtHistAnim.gif", fps = 15)

ftimeIdxs = filter(idx -> length(yInTime[idx]) > 1, timeIdxs)
speedYHistAnim = @animate for tidx in ftimeIdxs
    histogram2d(yInTime[tidx]*1e3, speedsInTime[tidx]/1e7,
                xlabel = "Position along channel [mm]",
                ylabel = "Speed [10⁷ m/s]",
                title = "t ≈ $(round(particles[1].times[tidx]*1e12, digits = 1))ps")
end
mp4(speedYHistAnim, "/tmp/speedYHist-1e4.mp4")
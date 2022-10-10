using DataFrames
using CSV
using CairoMakie
using Statistics
using StatsBase
using Printf

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
NstdErr(uₚ, f) = λ*uₚ/c/f/h

rationalFloat(x::Float64) = Float64(rationalize(x))
normalizedNumber(x::Float64) = modf(x)[1] == 0. ? Int(x) : x

function roundToError(value::Float64, error::Float64)
    roundedError = round(error, sigdigits = 1)
    multilpier = exp10(floor(log10(roundedError)))
    roundedValue = rationalFloat(round(value / multilpier) * multilpier)

    normalizedNumber(roundedValue), normalizedNumber(roundedError)
end

function SmoothedZscoreAlgo(y, lag, threshold, influence)
    # Julia implimentation of http://stackoverflow.com/a/22640362/6029703
    n = length(y)
    signals = zeros(n) # init signal results
    filteredY = copy(y) # init filtered series
    avgFilter = zeros(n) # init average filter
    stdFilter = zeros(n) # init std filter
    avgFilter[lag - 1] = mean(y[1:lag]) # init first value
    stdFilter[lag - 1] = std(y[1:lag]) # init first value
    
    for i in range(lag, stop=n-1)
        if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1]
            if y[i] > avgFilter[i-1]
                signals[i] += 1 # postive signal
            else
                signals[i] += -1 # negative signal
            end
            # Make influence lower
            filteredY[i] = influence*y[i] + (1-influence)*filteredY[i-1]
        else
            signals[i] = 0
            filteredY[i] = y[i]
        end
        avgFilter[i] = mean(filteredY[i-lag+1:i])
        stdFilter[i] = std(filteredY[i-lag+1:i])
    end
    return (signals = signals, avgFilter = avgFilter, stdFilter = stdFilter)
end

data = DataFrame(CSV.File("data/DATA22.CSV", header = [:P, :t], skipto = 5))
bkg = DataFrame(CSV.File("data/DATA21.CSV", header = [:P, :t], skipto = 5))
outData = DataFrame(CSV.File("freq", header=[:hex, :f]))

# try to find step edges
pdiff = abs.(diff(data.P))
res = SmoothedZscoreAlgo(data.P, 30, 4, 10)

#= visualy check the edges
plt = Figure()
ax = Axis(plt[1, 1])
lines!(ax, data.P ./ maximum(data.P))
lines!(ax, abs.(res.signals))
plt
=#


sigInd = findall(s -> abs(s) == 1.0, res.signals)
stepEdges = sort([0, sigInd[1:end-1]..., sigInd[6] + 67, sigInd[7] + 54, sigInd[end-1] + 55, length(data.P)])

# split data
steps = AbstractArray[]
for i in 2:length(stepEdges)
    push!(steps, data.P[stepEdges[i-1]+1:stepEdges[i]])
end

#clean data
meanBkg = mean(bkg.P)
cleanSteps = AbstractArray[]
for (i, step) in enumerate(steps)
    while true
        m = mean(step)
        sd = std(step, mean = m)

        garbage = findall(x -> abs(m - x) >= 3*sd, step)
        isempty(garbage) && break

        println("Removing $garbage in step $i")
        step = deleteat!(step, garbage)
    end

    push!(cleanSteps, step .- meanBkg)
end

outData[!, :P] = mean.(cleanSteps)
outData[!, :stderr] = std.(cleanSteps)
outData[!, :N] = N.(outData.P, outData.f .* 1e6)
outData[!, :Nerr] = NstdErr.(outData.stderr, outData.f .* 1e6)
outData[!, :NW] = Weights(1 ./ outData.Nerr .^2)

meanN = mean(outData.N, outData.NW)
meanNErr = std(outData.N, outData.NW, mean = meanN)

#= visualy check the means
plt = Figure()
ax = Axis(plt[1, 1])
lines!(ax, data.t, data.P)
hlines!.(ax, outData.P)
plt
=#

# P vs F plot
pltPF = Figure()
axPF = Axis(pltPF[1, 1], xlabel = "f [MHz]", ylabel = "P [μW]", xscale = log10)
scatterlines!(axPF, outData.f, outData.P .* 1e6)
errorbars!(axPF, outData.f, outData.P, outData.stderr, whiskerwidth = 10)

# N vs F plot
pltNF = Figure()
axNF = Axis(pltNF[1, 1], xlabel = "f [MHz]", ylabel = "Number of photons", limits=(nothing, nothing, 0, nothing), xscale = log10)
scatterlines!(axNF, outData.f, outData.N, color = :blue)
errorbars!(axNF, outData.f, outData.N, outData.Nerr, whiskerwidth = 10, color = :blue)

meanNRounded, roundedErr = roundToError(meanN, meanNErr)
hlines!(axNF, meanNRounded, color = :orange, label = "Mean: $(@sprintf("%.0e", meanNRounded)) ± $(@sprintf("%.0e", roundedErr))")
axislegend(axNF)
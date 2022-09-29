import JSON
using Trapz
using DelimitedFiles

struct Parameters
    Vs::Float64     # V
    C::Float64      # F
    R::Float64      # Ω
    Lmcp::Float64   # m
    k::Float64      # m⁻¹
    G::Float64
    Qs::Float64     # C
    L::Float64
    vs::Float64     # V
    RC::Float64

    Parameters(dict::Dict{Symbol, Any}) = begin
        args = [Float64(dict[fn]) for fn in fieldnames(Parameters)[1:end-5]]

        G = :G in keys(dict) ? Float64(dict[:G]) : log(args[1]/Float64(dict[:Vm]))
        Qs = args[1]*args[2]
        L = args[5]*args[4]
        vs = args[1]/L
        RC = args[2]*args[3]

        new(args..., G, Qs, L, vs, RC)
    end
end

# function to parse times and xcoords
function parseDomain(obj::Dict{Symbol, Any})
    if obj[:type] == "list"
        return Float64.(obj[:value])
    elseif  obj[:type] == "range"
        return eval(Meta.parse("range($(Float64(obj[:min])), $(Float64(obj[:max])), $(:len in keys(obj) ? "length=$(Int(obj[:len]))" : "step=$(Float64(obj[:step]))"))"))
    end

    throw(DomainError(obj[:type], "Unrecognized type"))
end

# function to parse i0 and psi0
function parseInput(obj::Dict{Symbol, Any}, domain)
    type = obj[:type]

    if type == "constant"
        return fill(Float64(obj[:value]), length(domain))
    elseif  type == "list"
        return Float64.(obj[:value])
    elseif type == "function"
        expr = Meta.parse(obj[:expr])
        f = eval(expr)
        return Base.invokelatest.(f, domain)
    end

    throw(DomainError(type, "Unrecognized type"))
end

# redefine log function to not throw error when argument is negative, but return NaN
function ntlog(arg)
    if arg < 0
        @warn "Negative logarithm argument: $arg"
        return NaN
    end

    log(arg)
end

function calculateGain(par::Parameters, i₀, ψ₀, times, xcoords)
    tsize = length(times)
    xsize = length(xcoords)
    tind = 1:tsize
    xind = 1:xsize

    ϵs = [range(0.1, 1, length = 10); fill(1, 20)]
    N = length(ϵs)
    ϵind = 1:N

    gs = Array{Float64}(undef, xsize, tsize, N+1)
    QQss = Array{Float64}(undef, xsize, tsize, N)
    QwQss = Array{Float64}(undef, tsize, N)
    ψs = Array{Float64}(undef, xsize, tsize, N)

    gs[:, :, 1] = [exp(par.G*x) for x in xcoords, t in tind] # g₀

    # calculate the gain
    for i in ϵind
        g = gs[:, :, i]

        QQs = [trapz(times[1:t], exp.(times[1:t] ./par.RC) .* i₀[1:t] .* g[x, 1:t])/par.Qs for x in xind, t in tind]
        QQss[:, :, i] = QQs

        QwQs = [trapz(xcoords, QQs[:, t])/par.L for t in tind]
        QwQss[:, i] = QwQs

        ψ = [(ψ₀[x] + (QwQs[t] - QQs[x, t])*ϵs[i])*exp(-times[t]/par.RC) for x in xind, t in tind]
        ψs[:, :, i] = ψ

        new_g = [exp(par.G*xcoords[x] + trapz(xcoords[1:x], ntlog.(1 .+ ψ[1:x, t]))) for x in xind, t in tind]
        gs[:, :, i+1] = new_g
    end

    gs, QQss, QwQss, ψs
end

function dumpTensor(tensor, steps, out::String; array = false)
    mkpath(out)

    for n in steps
        path = joinpath(out, "$n.csv")
        array ? writedlm(path, tensor[:, n]) : writedlm(path, tensor[:, :, n])
    end
end

function dumpTensor(tensor, out::String; array = false)
    dumpTensor(tensor, 1:size(tensor)[end], out; array = array)
end

# MAIN
conf = JSON.parsefile("sim-param.json", dicttype = Dict{Symbol, Any})

params = Parameters(conf)
xcoords = parseDomain(conf[:xcoords])
times = parseDomain(conf[:times])

i₀ = parseInput(conf[:i0], times)
ψ₀ = parseInput(conf[:psi0], xcoords)

gs, QQss, QwQss, ψs = calculateGain(params, i₀, ψ₀, times, xcoords)

# output
tesnsors = ["gs", "QQss", "ψs"]
translate = Dict("i0" => i₀, "psi0" => ψ₀)

outputDir = String(conf[:outputDir])
mkpath(outputDir)

writedlm(joinpath(outputDir, "times.csv"), times)
writedlm(joinpath(outputDir, "xcoords.csv"), xcoords)
write(joinpath(outputDir, "parameters.json"), JSON.json(params, 2))

if conf[:output] == "all"
    writedlm(joinpath(outputDir, "i0.csv"), i₀)
    writedlm(joinpath(outputDir, "psi0.csv"), ψ₀)

    for tensor in tesnsors
        dumpTensor(eval(Meta.parse(tensor)), joinpath(outputDir, tensor))
    end

    dumpTensor(QwQss, joinpath(outputDir, "QwQss"); array = true)

else
    for var in String.(conf[:output])
        parts = split(var, "&")

        if parts[1] in tesnsors || parts[1] == "QwQss"
            array = parts[1] == "QwQss"

            if length(parts) > 1
                dumpTensor(eval(Meta.parse(parts[1])), eval(Meta.parse(parts[2])), joinpath(outputDir, parts[1]); array = array)
            else
                dumpTensor(eval(Meta.parse(parts[1])), joinpath(outputDir, parts[1]); array = array)
            end

        else
            writedlm(joinpath(outputDir, "$(parts[1]).csv"), translate[parts[1]])
        end
    end
end

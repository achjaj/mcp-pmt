mutable struct Particle
    pidx::Int
    times::Vector{Float64}
    features::Dict{AbstractString, Vector{Float64}}
    units::Dict{AbstractString, AbstractString}
end

function parseExport(path::AbstractString, delim = '\t')
    pattern = r"\s([^\s]*)"
    particles = Particle[]

    io = open(path, "r")

    # skip header
    for i = 1:8
        readline(io)
    end

    # parse indicies
    indicies = parse.(Int, split(readline(io), delim))
    for i in indicies
        push!(particles, Particle(i, Float64[], Dict{AbstractString, Vector{Float64}}(), Dict{AbstractString, AbstractString}()))
    end

    feature = nothing
    while !eof(io)
        line = readline(io)

        if line[1] == '%'
            matches = collect(eachmatch(pattern, line))

            feature = matches[1].captures[1]
            unit = matches[2].captures[1]
            time = parse(Float64, split(matches[4].captures[1], "=")[end])

            for particle in particles
                (length(particle.times) == 0 || particle.times[end] != time) && push!(particle.times, time)
                
                if !haskey(particle.features, feature)
                    particle.features[feature] = Float64[]
                    particle.units[feature] = unit
                end
            end

        else
            values = parse.(Float64, split(line, delim))
            
            for (i, value) in enumerate(values)
                push!(particles[i].features[feature], value)
            end
        end
    end

    close(io)
    return particles
end
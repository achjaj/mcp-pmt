using PlotlyJS

function skiplines(io::IO, n::Int64)
    for i in 1:n
        readline(io)
    end
end

function readQySection(io::IO, ln::Int64)
    line = readline(io)
    ch = line[1]
    ch != '%' && throw(ErrorException("Line $(ln + 1): Expected start of section"))

    time = -1.0
    coordinates = Float64[]

    if occursin("qy", line)
        time = parse(Float64, split(split(line, "@")[2], "=")[2])

        line = readline(io)
        coordinates = parse.(Float64, split(line, "\t"))
    else
        readline(io)
    end

    return time, coordinates
end


function parseFile(path::String)
    times = Float64[]
    #x = Vector{Vector{Float64}}()
    y = Vector{Vector{Float64}}()
    #z = Vector{Vector{Float64}}()


    io = open(path, "r")

    skiplines(io, 9)
    ln = 9

    while !eof(io)
        time, coords = readQySection(io, ln)

        if time > -1.0
            push!(times, time)
            push!(y, coords)
        end

        ln += 2
    end

    close(io)

    return times, y
end

times, y = parseFile("position.csv")

by_particle = Vector{Vector{Float64}}()

for i in 1:length(y)
    push!(by_particle, getindex.(y, i))
end

diffs = diff.(by_particle)
filtered_diffs = filter.(x -> !isnan(x), diffs)
max_diffs = maximum.(filtered_diffs, init = -1)
max_diff = maximum(max_diffs)

L = 4.2e-4
counts = Int64[]
maxY = Float64[]

for coords in y
    push!(counts, length(filter(x -> (L - max_diff) <= x <= (L + max_diff), coords)))
    
    filtered = filter(x -> !isnan(x), coords)
    if isempty(filtered)
        push!(maxY, NaN)
    else 
        push!(maxY, maximum(filtered))
    end
end

plot(scatter(x = times, y = counts))
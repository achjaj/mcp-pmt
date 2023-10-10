using PlotlyJS

counts = Int64[]

open("counts.csv", "r") do io
    for i in 1:8
        readline(io)
    end
    
    global times
    times = split(readline(io), "\t")[4:end]
    times = parse.(Float64, getindex.(split.(times, "t="), 2))

    line = split(readline(io), "\t")[4:end]

    global counts
    counts = parse.(Int64, line)
end

plot(scatter(x = times, y = counts, type = "line"), Layout(yaxis_type = "log", xaxis_title = "Time [s]", yaxis_title = "Number of electrons"))
using GLMakie
using DelimitedFiles

import CairoMakie

const cmdl = "\$> "

# define commands
Base.@kwdef struct Plot 
    help = "Plot data"
end

Base.@kwdef struct Export 
    help = "Export plot"
end

Base.@kwdef struct ListSteps
    help = "List available step data"
end

Base.@kwdef struct Show
    help = "Print file content"
end

Base.@kwdef struct Help
    help = "Print this help"
end

cmds = ["Plot", "Export", "ListSteps", "Show", "Help"]

run_cmd(x::ListSteps, dataName) = (i -> println(i), readdir(dataName))
run_cmd(x::Show, fileName) = map(l -> println(l), readlines(fileName))
function run_cmd(x::Help, args...)
    println("Available commands:")
    for cmd in cmds
        println("\t$cmd - $(eval(Meta.parse("$cmd().help")))")
    end
end

function run_cmd(x::Plot, dataName, )

function initCheck()
    lookFor = ["xcoords.csv", "times.csv"]

    for f in lookFor
        isfile(f) || throw(ErrorException("Could not find '$f' in TLM output directory!"))
    end
end

# change working directory to supplied TLM output directory
cd(ARGS[1])
#cd("output/")

initCheck()
xcoords = readdlm("xcoords.csv")
times = readdlm("times.csv")

# main loop
while true
    println("\nAvailable data:\n$(readdir())")
    println("(Type Help for help)")
    print(cmdl)

    splitLine = split(readline(), " ")
    cmdLine = Array{Any}(undef, length(splitLine))
    
    cmdLine[1] = eval(Meta.parse("$(splitLine[1])()"))
    cmdLine[2:end] = splitLine[2:end]

    run_cmd(cmdLine...)
end
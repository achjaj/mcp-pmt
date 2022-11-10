using Pkg

# check required packages
println("Checking for dependencies")
pkgs = ["Distributions", "GLMakie", "GeometryBasics"]

# get list of installed packages
buffer = IOBuffer()
Pkg.status(io = buffer)
buffer.ptr = 1
installed = [line[end-1] for line in split.(readlines(buffer)[2:end-1], " ")]

# install missing packages
mis = filter(pkg -> !(pkg in installed), pkgs)
isempty(mis) || Pkg.add.(mis)

println("Continuing with the simulation")
include("kinematics.jl")
include("drawing.jl")
include("see.jl")

function oneArc(T, y₀)
    println("\t\tT: $T")
    v0 = rndVec(v(T*q)) # the energy is in eV, but we want the velocity in SI
    t = tc(v0[[1, 3]])

    return v0, t, r(v0, t)[2] + y₀
end

function oneElectron(queue, T₀ = 5, y₀ = 0)
    println("\tStep 1")
    steps = [oneArc(T₀, y₀)]

    sn = 2;
    nT₀ = T₀
    while steps[end][end] < L
        println("\tStep $sn")
        s = steps[end]

        vc = v̅(s[1:2]...)
        θc = asin(vc[2]/norm(vc))

        S = length(steps) == 1 ? s[end] : s[end] - steps[end-1][end]
        Tc =  V*S/L + nT₀ #1/2*m*norm(vc)^2 / q # in eV

        println("\t\tTc: $Tc")

        seeE = see(Tc, θc)
        println("\t\tseeE: $seeE")

        #length(seeE) == 0 && break

        #=for energy in seeE[2:end]
            push!(queue, (energy, s[end]))
        end=#

        nT₀ = rand(Rayleigh(5))#length(seeE) == 1 ? seeE[1] : Tc - sum(seeE)

        push!(steps, oneArc(nT₀, s[end]))

        sn += 1
    end

    steps
end

# run the simulation
function run(;V_in = 2e3, L_in = 1.5e-3, R_in = 1e-5, q_in = 1.6e-19, m_in = 9.1e-31)
    # MCP-PMT parameters
    global V
    global L
    global R
    global q
    global m
    global a

    V = V_in
    L = L_in
    R = R_in
    q = q_in
    m = m_in

    # e⁻ acceleration
    a = q*V/L/m


    queue = [(5.0, 0.0)] # queue of electrons
    result = []

    en = 1
    #while !isempty(queue)
        println("Working on electron $en")
        sv = popfirst!(queue)
        push!(result, oneElectron(queue, sv...))

        en += 1
    #end

    return result
end

result = run()
steps = result[1]
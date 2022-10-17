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

println("Continuing with the initialization")#
using Distributions
using GLMakie
using GeometryBasics
using LinearAlgebra

# MCP-PMT params
V = 2400    # V
L = 1.5e-3  # m
R = 1e-6    # m
E = V/L     # V/m
q = 1.6e-19 # C
m = 9.1e-31 # kg
k = 1/2

# e⁻ acceleration
a = q*E/m

# initialize distributions for generating random values
ϕDist = Cosine(π, π)        # ϕ ∈ <0, 2π>
θDist = Cosine(π/4,π/4)     # θ ∈ <0, π/2>
v₀Dist = Rayleigh(2e5)     # distribution of magnitude of initial velocity

# generate random vector with magnitude mag; default magitude is 1
function rndVec(mag = 1.0)
    # generate random angles
    ϕ = rand(ϕDist)
    θ = rand(θDist)

    # return the vector in cartesian coordinates and the angles
    return mag*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)], ϕ, θ
end

# calculates the path of an e⁻
r(v₀, t) = v₀*t + [0, a*t^2, 0]

# calculates the kinetic energy of e⁻
T(v₀, t) = m*(norm(v₀*t + [0, a*t, 0])^2)/2

# calculates the time of collision
tc(v) = 2R*v[3]/(v[1]^2 + v[3]^2)

# creates plot with a cylinder, i. e. our channel
function preparePlot()
    #creates empty plot
    plt = Figure()
    ax = Axis3(plt[1, 1], aspect = (1, 3, 1))

    # draw the cylinder in red color and transparency
    mesh!(ax, Cylinder{3, Float32}(Point3f0(0.0, 0.0, R), Point3f0(0.0, L, R), R), color = (:red, 0.3), transparency=true)
    
    return ax, plt
end


function drawSteps(steps, res = 100)
    ax, plt = preparePlot()

    l = 0
    d_glob = zeros(3)
    for s in steps
        d = r(s[1:2]...)
        
        l += 2*R*asin(norm(d[[1, 3]])/R/2)
        if l > 2π*R
            l -= R
        end

        β = l/R

        d_glob_new = [R*cos(β), d_glob[2] + d[2], R + R*sin(β)]
        v = ([d_glob_new[1], d_glob_new[2] - a*s[2]^2, d_glob_new[3]] - d_glob) / s[2]

        times = range(0, s[2], length = res)
        lines!(ax, hcat([r(v, t) + d_glob for t in times]...), color = :blue)

        d_glob = d_glob_new
    end

    plt, ax
end

# calculates one step
function step()
    magv₀ = rand(v₀Dist)
    v₀, ϕ, θ = rndVec(magv₀)
    t = tc(v₀)                 # time of collision

    Tc = T(v₀, t)              # collision energy
    δ = k*Tc

   return v₀, t, δ, r(v₀, t)[2], ϕ, θ
end

# run the simulation and draw plot of the paths
function run()
    steps = Tuple[]
    n = 0
    while n < 1 || sum([step[4] for step in steps]) <= L
        n += 1
        newStep = step()
        push!(steps, newStep)
        println("Step $n:\n\tS = $(sum([s[4] for s in steps]))")
    end

    plt, ax = drawSteps(steps)
    display(plt)

    return steps, plt
end

run()
using Pkg

# check required packages
println("Checking for dependencies")
pkgs = ["Distributions", "GLMakie", "GeometryBasics", "LinearAlgebra"]

# get list of installed packages
buffer = IOBuffer()
Pkg.status(io = buffer)
buffer.ptr = 1
installed = [line[end-1] for line in split.(readlines(buffer)[2:end-1], " ")]

# install missing packages
mis = filter(pkg -> !(pkg in installed), pkgs)
isempty(mis) || Pkg.add.(mis)

println("Continuing with the initialization")
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

# initialize distributions for generating random values
ϕDist = Cosine(π, π)        # ϕ ∈ <0, 2π>
θDist = Cosine(π/4,π/4)     # θ ∈ <0, π/2>
T₀Dist = Rayleigh(1e-20)

# generate random vector with magnitude mag; default magitude is 1
function rndVec(mag = 1.0)
    # generate random angles
    ϕ = rand(ϕDist)
    θ = rand(θDist)

    # return the vector in cartesian coordinates and the angles
    return mag*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)], ϕ, θ
end

# test the generation of random vectors; creates 3D plot of vector endpoints
function rndVectorTest(n::Int = 10000)
    vcs = map(x -> getindex(rndVec(), 1), 1:n)

    vcx = map(v -> getindex(v, 1), vcs)
    vcy = map(v -> getindex(v, 2), vcs)
    vcz = map(v -> getindex(v, 3), vcs)

    meshscatter(vcx, vcy, vcz, color = :red, markersize = 1e-2)
end

# calculates rotatin matrix
function rotMat(xz)
    β = -2*asin(norm(xz)/(2*R))

    return [cos(β) 0 sin(β)
                0  1 0
            -sin(β) 0 cos(β)]
end

# calculates the path of an e⁻
xz(v₀, t) = v₀*t
y(v₀, t) = v₀*t + q*E*(t^2)/m

# calculates the speed of e⁻ in y direction
vy(v₀, t) = v₀ + q*E*t/m

# calculates the kinetic energy of e⁻
T(v₀, t) = m*(v₀[1]^2 + vy(v₀[2], t)^2 + v₀[3]^2)/2

# calculates the time of collision
tc(vx, vz) = 2R*vz/(vx^2 + vz^2)

# creates plot with a cylinder, i. e. our channel
function preparePlot()
    #creates empty plot
    plt = Figure()
    ax = Axis3(plt[1, 1], aspect = (1, 3, 1))

    # draw the cylinder in red color and transparency
    mesh!(ax, Cylinder{3, Float32}(Point3f0(0.0, 0.0, R), Point3f0(0.0, L, R), R), color = (:red, 0.3), transparency=true)

    return ax, plt
end

# draw one step
function drawStep(ax, s, rot, r₀; res = 100)
    times = range(0, s[3], length = res)
    v₀ = s[2]

    # calculates the path and converts it to global coordinate system
    vec = hcat([rot*[xz(v₀[1], t), y(v₀[2], t), xz(v₀[3], t)] + r₀ for t in times]...)

    # draw the path
    lines!(ax, vec[1, :], vec[2, :], vec[3, :], color = :blue)

    # highlight start and collision position: green dot is the start position, red is the collision position
    scatter!(ax, vec[1, [1, end]], vec[2, [1, end]], vec[3, [1, end]], color = [:green, :red])
end

# draw all the steps
function drawSteps(ax, steps)
    globr₀ = zeros(3)

    for (i, s) in enumerate(steps)
        println("Drawing step $i; r₀: $globr₀")
        rot = rotMat(globr₀[[1, 3]])
        drawStep(ax, s, rot, globr₀)

        v₀ = s[2]
        t = s[3]
        S = s[4]
        
        globr₀ = rot*([xz(v₀[1], t), S, xz(v₀[3], t)] * globr₀)
    end
end

# calculates one step
function step()
    T₀ = rand(T₀Dist)           # initial energy
    magv₀ = sqrt(2*T₀/m)
    v₀, ϕ, θ = rndVec(magv₀)
    t = tc(v₀[1], v₀[3])        # time of collision

    S = y(v₀[2], t)
    Tc = T(v₀, t)              # collision energy
    δ = k*Tc

   return T₀, v₀, t, S, Tc, δ, ϕ, θ
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

    ax, plt = preparePlot()
    drawSteps(ax, steps)
    display(plt)

    return steps, plt
end

run()
#import GLMakie as GL 
using GeometryBasics
using LinearAlgebra
import Plots
Plots.plotlyjs()

include("parseExport.jl")

channelLen = 4.2e-4 # m
U = 2e3
q = 1.602176634e-19
m = 9.1093837015e-31
ltd = 42
R = channelLen / 42 / 2
boundary = channelLen
fieldAngle = deg2rad(0) #rad

a = [1 0 0; 0 cos(fieldAngle) -sin(fieldAngle); 0 sin(fieldAngle) cos(fieldAngle)] * [0, q*U/channelLen/m, 0]

lastYidx(y) = findmax(value -> isnan(value) ? -Inf : value, y)[end]
calcTheta(vy, speed) = acos(vy/speed) # calculates the angle between the velocity vector and the normal to the channel opening face; the normal is vector [0, 1, 0]

function intersects(pos, velocity)
    vy = velocity[2]
    ay = a[2]
    y = pos[2]

    t = (-vy + sqrt(vy^2 + 2*ay*(channelLen-y)))/ay

    p = pos + velocity*t + 0.5*a*t^2

    p[1]^2 + p[3]^2 <= R^2
end

getVelocity(particle, i) = [particle.features["cpt.vx"][i], particle.features["cpt.vy"][i], particle.features["cpt.vz"][i]]
getPos(particle, i) = [particle.features["qx"][i], particle.features["qy"][i], particle.features["qz"][i]]

particles = parseExport("pivi/particle-pivi-1e4-biased_field_15deg-all.txt")
#particles = parseExport("pivi/particle-pivi-1e4-all.txt")

x = Float64[]
y = Float64[]
z = Float64[]

u = Float64[]
v = Float64[]
w = Float64[]

lengths = Float64[]
thetas = Float64[]

phis = Float64[]

for particle in particles
    i = lastYidx(particle.features["qy"])
    velocity = getVelocity(particle, i)
    position = getPos(particle, i)

    if intersects(position, velocity)
        push!(x, position[1])
        push!(y, position[2])
        push!(z, position[3])

        push!(u, velocity[1])
        push!(v, velocity[2])
        push!(w, velocity[3])

        speed = norm(velocity)
        push!(lengths, speed)
        push!(thetas, calcTheta(velocity[2], speed))

        push!(phis, atan(position[3], position[1]))
    end
end

energies = 0.5 .* m .* lengths.^2 ./ q # in eV
thetas = rad2deg.(thetas)
filter!(x -> x <= 15, thetas)

#=fig = GL.Figure()

axis, arrows_obj = GL.arrows(fig[1, 1], x, y, z, u, v, w, 
        fxaa = true,
        axis=(type=Axis3,),
        linewidth = 5e-7,
        arrowsize = GL.Vec3(2e-6, 2e-6, 3e-6),
        lengthscale = 1e-12,
        color = lengths
    )=#

#mesh!(axis, Cylinder{3, Float32}(Point3f0(0.0, channelLen - channelLen/20, R), Point3f0(0.0, channelLen, R), R), color = (:red, 0.3), transparency=true)

#GL.Colorbar(fig[1, 2], limits = extrema(lengths), label = "Speed [m/s]")

#=noFrames = 3*120
GL.record(fig, "pivi/outVel.mp4", 1:noFrames) do frame
    axis.azimuth[] = 2π/noFrames * frame
end=#

#speedHist = Plots.stephist(lengths/1e7, xlabel = "Exit speed [10⁷ m/s]", ylabel = "Count", legend = false)
engHist = Plots.stephist(energies ./ 1000, xlabel = "Exit energy [keV]", ylabel = "Count", legend = false)
thetaHist = Plots.stephist(thetas, xlabel = "θ [degree]", ylabel = "Count", legend = false)
#=speedThetaHist = Plots.histogram2d(thetas, lengths/1e7, xlabel = "θ [rad]", ylabel = "Exit speed [10⁷ m/s]")
yHist = Plots.stephist(y .*1e3, xlabel = "Position along channel [mm]", ylabel = "Count", legend = false)
yPhiHist = Plots.histogram2d(y .*1e3, phis, xlabel = "Position along channel [mm]", ylabel = "Φ [rad]")
speedPosHist = Plots.histogram2d(y *1e3, lengths / 1e7, xlabel = "Position along channel [mm]", ylabel = "Exit speed [10⁷, m/s]")=#
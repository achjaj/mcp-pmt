using LinearAlgebra
using Distributions

# MCP-PMT params
V = 2e3     # V
L = 1.5e-3  # m
R = 1e-5    # m
q = 1.6e-19 # C
m = 9.1e-31 # kg

# e⁻ acceleration
a = q*V/L/m

# initialize distributions for generating random values
ϕDist = Uniform(0, π)      # ϕ ∈ [0, π]
θDist = Cosine(π/4,π/4)     # θ ∈ [0, π/2]

# generate random vector with magnitude mag; default magitude is 1
function rndVec(mag = 1.0)
    # generate random angles
    ϕ = rand(ϕDist)
    θ = rand(θDist)

    # return the vector in cartesian coordinates and the angles
    return mag*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

# calculates the path of an e⁻
r(v₀, t) = v₀*t + [0, (a*t^2)/2, 0]

# calculates the velocity at time t
v̅(v₀, t) = v₀ + [0, a*t, 0]

# calculates the time of collision
#tc(v_xz) = 2R*v_xz[2]/norm(v_xz)^2
#=function tc(v_xz)
    k = v_xz[2]/v_xz[1]

    2*sqrt(2)*R/v_xz[1]*k/(1+k^2)^(3/2)
end=#

function tc(v_xz)
    k = norm(v_xz)^2
    
    2R*v_xz[2]/k
end

# calculates speed from kinetic energy (in J)
v(T) = sqrt(2*T/m)
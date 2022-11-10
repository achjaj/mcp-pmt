include("kinematics.jl")
include("see.jl")
include("drawing.jl")

function moveAndRotate(vector, to)
    d̄ = to[[1, 3]]
    d = norm(d̄)

    ratio = d/R/2
    sgn = to[1] > 0 ? -1 : 1

    β = 2*asin(ratio)*sgn

    mat = rotMat(β)

    mat*vector + to
end

collided(r̄) = r̄[1]^2 + (r̄[3] - R)^2 >= R^2

function run(;V_in = 2e3, L_in = 1.5e-3, R_in = 1e-5, q_in = 1.6e-19, m_in = 9.1e-31, dt_in = 1e-15)
    # MCP-PMT parameters
    global V
    global L
    global R
    global q
    global m
    global a
    global Ē
    global dt

    V = V_in
    L = L_in
    R = R_in
    q = q_in
    m = m_in

    # e⁻ acceleration
    a = q*V/L/m

    Ē = [0, V/L, 0]
    dt = dt_in

    positions = [zeros(3)]
    velocities = [rndVec(v(5*q))]

    while positions[end][2] < L
        F̄ = q*Ē
        dv̄ = F̄/m*dt

        v̄ = velocities[end] + dv̄
        dr̄ = v̄*dt

        r̄ = positions[end] + dr̄

        if collided(r̄)
            speed = norm(v̄)
            Tc = 0.5*m*speed^2 / q
            θc = asin(v̄[2]/speed)

            T₀ = rand(Rayleigh(5)) #Tc - sum(see(Tc, θc))

            v̄ = moveAndRotate(rndVec(v(T₀*q)), r̄)
        end

        push!(positions, r̄)
        push!(velocities, v̄)
    end

    return positions, velocities
end

positions, velocities = run()
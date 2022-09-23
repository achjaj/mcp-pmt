using Trapz

#= matrices: 
    ---> t
    |
    |
    v x
=#

struct Params
    Vs::Number          # V
    R::Number           # Ω
    C::Number           # F
    Lmcp::Number        # m
    G::Number
    k::Number           # m^-1
end

# calculate gain in place x and time t using the iterative process
# parameters:
#   i₀ - input signal matrix: first column is time; second column is value i₀; the gain is calculated at the last time value
#   ψ₀ - ψ matrix at time 0: first column is coordinate x; second column is value; tha gain is calculated in the last x value
#   params - parameters of MCP-PMT
function findGain(i₀, ψ₀, params)
    Qs = params.Vs*params.C   # C
    Q₀ = trapz(i₀[:, 1], i₀[:, 2])
    L = params.k * params.Lmcp
    xs = collect(range(0, L, length = size(ψ₀)[1]))
    x = ψ₀[:, 1][end]


    gs = []
    push!(gs, exp.(params.G .* xs))
    for ϵ in fill(1, 30)
        Q = gs[end] .* Q₀
        Qw0 = (size(Q) == () ? trapz([xs[end]], [Q / Qs]) : trapz(xs, Q ./ Qs)) ./ L - Q₀/Qs
        ψ = ψ₀[:, 2] .+ ϵ*((Qw0/Qs) .+ (Q ./ Qs)) .+ (Q₀/Qs)
        g = exp.(params.G*x .+ trapz(ψ₀[:, 1], log.(1 .+ ψ)))
        
        push!(gs, g)
    end

    gs
end

params = Params(2400, 6.09e8, 2.86e-13, 1.5e-3, 0.78, 1.16e4)

i₀ = [range(0, 10e-9, length = 10) fill(1.75e-5, 10)]
ψ₀ = [range(0, 1e-3, length = 10) fill(0, 10)]

gs = findGain(i₀, ψ₀, params)
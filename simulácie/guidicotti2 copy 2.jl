using Trapz
using Plots

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

params = Params(2400, 6.09e8, 2.86e-13, 1.5e-3, 0.78, 1.16e4)

Qs = params.C * params.Vs
L = params.k * params.Lmcp

function sym(xi, Q₀, g₀, ψ₀, xgrid)
    
    gs = [g₀]
    for i in [range(0.1, 1, length = 10); ones(20)]
        g = gs[end]
  
        Q = (Q₀ .*g)/Qs
        Qw = Q - Q₀/Qs
        ψ = ψ₀ .+ (Q₀/Qs) .+ Qw .- Q
        g = exp.(params.G*xgrid[xi] + trapz(xgrid[1:xi], log.(1 .+ ψ[1:xi])))
        
        push!(gs, g)
    end

    gs
end

xgrid = range(0, L, length = 100)
xi = 50
g₀ = exp(params.G*xgrid[xi])
ψ₀ = zeros(length(xgrid))
Q₀ = (1.75e-5)*(20e-9)

gs = sym(xi, Q₀, g₀, ψ₀, xgrid)
println(gs)

#scatter(grid[:, 2] ./ params.k, gs[end], yscale = :log10, legend = false)

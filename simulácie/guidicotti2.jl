using Trapz
using Plots

#= matrices: 
    ---> t
    |
    |
    v x
=#

# prepare the inpot signal function
coeffs = readdlm("input-signal-poly.csv", comments = true)
i₀F = Polynomial(coeffs[:, 1])

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

mesh = hcat(range(0, 30e-9, length = 1000), range(0, L, length = 1000)) # time and position discretization
times = range(0, 50e-9, length = 6)

i₀ = fill(0, 1000) #i₀F.(mesh[:, 1])
ψ₀ = fill(0, 1000)

Q₀ = trapz(mesh[:, 1], i₀)

gs = [exp.(params.G .* mesh[:, 2])]

for i in [range(0.1, 1, length = 10); ones(20)]
    Q = [trapz(mesh[:, 1], i₀ .* gs[end]) for j in mesh[:, 2]]
    Qw = trapz(mesh[:, 2], Q/Qs) /L - Q₀/Qs
    ψ = ψ₀ .+ Q₀/Qs .+ ((Qw) .- (Q ./ Qs))
    g = map(j -> exp(params.G*mesh[j, 2] + trapz(mesh[1:j, 2], log.(1 .+ ψ[1:j]))), 1:size(mesh)[1])

    push!(gs, g)
end

plot(mesh[:, 2] ./ params.k, gs[end], yscale = :log10, legend = false)

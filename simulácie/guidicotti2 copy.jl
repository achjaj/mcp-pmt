using Trapz
using Plots
using Polynomials
using DelimitedFiles

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

# prepare the inpot signal function
coeffs = readdlm("input-signal-poly.csv", comments = true)
i₀F = Polynomial(coeffs[:, 1])

params = Params(2400, 6.09e8, 2.86e-13, 1.5e-3, 0.78, 1.16e4)

Qs = params.C * params.Vs
L = params.k * params.Lmcp

mesh = hcat(range(1e-9, 50e-9, length = 1000), range(0.1, L, length = 1000)) # time and position discretization
mesh_times_ind = Int.(round.(range(1, size(mesh)[1], length = 6)))
times_ind = 1:6

i₀ = i₀F.(mesh[:, 1]) #fill(1.75e-5, 1000)
ψ₀ = fill(0, 1000)

Q₀ = [trapz(mesh[1:t, 1], i₀[1:t]) for t in mesh_times_ind]

g₀ = [exp.(params.G * x) for x in mesh[:, 2], t in mesh_times_ind]
gs = [g₀]

for i in [range(0.1, 1, length = 10); ones(20)]
    Q = gs[end] .* Q₀ #[trapz(mesh[:, 1], i₀ .* gs[end]) for j in mesh[:, 2]]
    Qw = trapz(mesh[:, 2], Q/Qs) /L - Q₀/Qs
    ψ = ψ₀ .+ Q₀/Qs .+ ((Qw) .+ (Q ./ Qs))*i
    g = map(j -> exp(params.G*mesh[j, 2] + trapz(mesh[1:j, 2], log.(1 .+ ψ[1:j]))), 1:size(mesh)[1])

    push!(gs, g)
end

#=for i in [range(0.1, 1, length = 10); ones(20)]
    g = gs[end]
    Q = [Q₀[t]*g[xi, t]/Qs for xi in range(1, size(g)[1]), t in times_ind]
    Qw = [trapz(mesh[:, 2], Q[:, t]) / L - Q₀[t]/Qs for t in times_ind]
    ψ = [ψ₀[x] + Q₀[t]/Qs + Qw[t] - Q[x, t] for x in 1:size(mesh)[1], t in times_ind]
    g = [exp(params.G*mesh[x, 2] + trapz(mesh[1:x, 2], log.(1 .+ ψ[1:x, t]))) for x in 1:size(mesh)[1], t in times_ind]

    push!(gs, g)
end=#

plot(mesh[:, 2] ./ params.k, gs[end][:, 2], yscale = :log10, legend = false)

using Trapz
using PlotlyJS

#= matrices: 
    ---> t
    |
    |
    v x
=#

# some detector properties
Vs = 2400               # V
R = 6.09e8              # Ω
C = 2.86e-13            # F
L = 1.5e-3              # m
G = 0.78
k = 1.16e4              # m^-1

Qs = Vs*C               # C

i₀ = fill(17e-8, 6)                                             # sampled input signal; square signal with value of 17e-8 A
st = range(0, 50e-9, length = 6)                                # time of samples in seconds

Q₀ = trapz(st, i₀)                                              # initial charge in C in time 50e-9 s

sx = range(0, k*L, length = 10)                                 # sampled x coordinate (for initial g)
g₀ = [exp(G*x) for x in sx]                                     # sampled initial gain as function iof space and time

ψ₀ = fill(0, length(sx))                                        # ψ at t = 0 as function of x

QpQs(g_val) =  g_val .* (Q₀ / Qs)
Qw0pQs(QpQs_val::Vector, x::Vector) = trapz(x, QpQs_val)/L - (Q₀/Qs)
Qw0pQs(QpQs_val::Number, x::Number) = QpQs_val - (Q₀/Qs)
ψ(QpQs_val, Qw0pQs_val, ϵ=1) = ϵ * (QpQs_val .+ Qw0pQs_val) .+ (Q₀/Qs) .+ ψ₀
g(ψ_val, x) = exp(G*x[end] + trapz(x, log.(1 .+ ψ_val)))

gs = [g₀]
xs = [sx]
ϵs = [range(0.1, 1, length = 10); fill(1, 20)]
xi = 1:1:length(sx)
for ϵ in ϵs
    println(gs[end])
    QpQs_val = QpQs(gs[end])
    Qw0pQs_val = Qw0pQs(QpQs_val, xs[end])
    ψ_val = ψ(QpQs_val, Qw0pQs_val)
    
    new_g = g(ψ_val, [k*(1e-3)])
    push!(gs, new_g)
    push!(xs, [k*(1e-3)])
end
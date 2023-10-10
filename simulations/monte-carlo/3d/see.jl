using Distributions
using SpecialFunctions

pn = [2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]
ϵn = [1.5, 1.75, 1, 3.75, 8.5, 11.5, 2.5, 3, 2.5, 3]

H(x) = Int(x > 0) # Heaviside step function

# here I rename some functions
Γ(x) = gamma(x)
P(a, x) = gamma_inc(a, x)[1]
invP(a, p) = gamma_inc_inv(a, p, 1 - p)
invβ(x, a, b) = beta_inc_inv(a, b, x)[1]
Π(arr) = *(arr...)

nSphere2Cart(R, ϕ) = [R*cos(ϕ[1]); [R*Π(sin.(ϕ[1:k-1])*cos(ϕ[k])) for k in 2:length(ϕ)]; R*Π(sin.(ϕ[1:length(ϕ)]))]

D(x, s = 1.54) = s*x/(s-1+x^s)
δ̅ts(θ; δ = 1.9, t1 = 0.7, t2 = 0.8) = δ*(1 + t1*(1 - cos(θ)^t2))
E̅ts(θ; E = 277, t3 = 0.7, t4 = 1) = E*(1 + t3*(1 - cos(θ)^t4)) 
δe0(E; P1e_inf = 0.02, P1e = 0.496, E̅e = 0, W = 60.86, p = 1) = P1e_inf + (P1e - P1e_inf)*exp(-((abs(E - E̅e)/W)^p)/p)
δr0(E; P1r_inf = 0.2, Er = 0.041, r = 0.104) = P1r_inf*(1 - exp(-(E/Er)^r)) 


f1e(E, E0, δe, σ = 2) = H(E)*H(E0 - E)*δe*(2*exp(-((E-E0)^2)/(2*σ^2)))/(sqrt(2π)*σ*erf(E0/(sqrt(2)*σ)))
f1r(E, E0, δr, q = 0.5) = H(E)*H(E0 - E)*δr*(q + 1)*(E^q)/E0^(q+1)


Fn(n, Pn, E0, ϵ, p) = Pn/((ϵ^p*Γ(p))^n * P(n*p, E0/ϵ))
fnts(E, F, p, ϵ) = H(E)*F*E^(p - 1)*exp(-E/ϵ)

function generateOneEnergy(δe::Number, δr::Number, P1::Number, E0; σ = 2, q = 0.5)
    δ1 = δe + δr + P1

    ae = δe/δ1
    ar = δr/δ1

    u = rand(Uniform()) # random number from [0, 1] with uniform distribution

    0 <= u < ae && return E0 - σ*abs(rand(Normal()))
    ae <= u < ae + ar && return E0*rand(Uniform())^(1/(1+q))
    ae + ar <= u < 1 && return ϵn[1] * invP(pn[1], rand(Uniform())*P(pn[1], E0/ϵn[1]))
end

function generateNEnergies(n::Integer, E0::Number)
    x0 = E0/ϵn[n]
    P0 = P(n*pn[n], x0)

    u = rand(Uniform(), n-1)
    μ = pn[1:n-1] .* (n .- 1:n-1) # μ = pₙ*(n - k), where k ∈ [1, n]
    ν = pn[1:n-1]
    θ = asin.(sqrt.(invβ.(u, μ, ν)))

    Y = sqrt(invP(n*pn[n], rand(Uniform()*P0)))
    y = nSphere2Cart(Y, θ)

    return ϵn[1:n] .* y.^2 
end

function samplePDF(pdf::AbstractArray)
    cdf = cumsum(pdf)
    a = rand() * maximum(cdf)

    a <= cdf[1] && return 1
    for i in 2:length(cdf)
        cdf[i - 1] < a <= cdf[i] && return i
    end
end

function getδs(T::Number, θ::Number; e1 = 0.26, e2 = 2, r1 = 0.26, r2 = 2)
    δe = δe0(T)*(1 + e1*(1 - cos(θ)^e2))
    δr = δr0(T)*(1 + r1*(1 - cos(θ)^r2))
    δts = δ̅ts(θ)*D(T/E̅ts(θ))

    δ̃ts = δts/(1 - δe - δr)

    return δe, δr, δ̃ts
end

function getPn(δe::Number, δr::Number, δ̃ts::Number, M::Integer)
    P̃n_ts = n -> pdf(Poisson(δ̃ts), n)

    Pn = (1 - δe - δr) .* P̃n_ts.(0:M)
    Pn[2] = Pn[2] + δe + δr

    return Pn
end

function see(T::Number, θ::Number; M = 10)
    δe, δr, δ̃ts = getδs(T, θ)

    Pn = getPn(δe, δr, δ̃ts, M)

    n = samplePDF(Pn) - 1
    
    n == 0 && return Float64[]
    
    n == 1 && return [abs(generateOneEnergy(δe, δr, Pn[2], T))]

    return generateNEnergies(n, T)
end
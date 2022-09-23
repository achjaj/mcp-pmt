using NLsolve

function simulate(i₀, Vₛ, R, G, L, k; a = 2.05, A₀ = 0.025, threshold = 1e-3)
    Iₛ = Vₛ/R
    B = i₀/Iₛ
    H = G + a
    vₛ = Vₛ/(k*L)

    C(A) = A + B
    D(A) = 2*G/(G+a) + C(A)
    E(A) = C(A) + 2

    function get_implicit_g(A, x)
        D_val = D(A)
        E_val = E(A)
        M = D_val - E_val

        function im_g!(F, g)
            F[1] = E_val*log(g[1]) + M*log((B*g[1] - D_val)/(B - D_val)) - D_val*H*x 
        end

        return im_g! 
    end

    nextA(ΦL, gL) = (G+a)/(2*vₛ*log(gL))*ΦL + 2*G*L/log(gL) - B - 2
    g(A, x) = nlsolve(get_implicit_g(A, x), [0.1]).zero[1]
    Φ(A, g_x, x) = vₛ*((C(A) - D(A))*x + E(A)/H*log(g_x))

    # firstly we need to find A using iterative process
    As = [A₀]
    ΦLs = []
    gLs = []

    while true
        println("Step $(length(As))")
        currA = As[end]

        print("\tCalculating gL: ")
        gL = g(currA, L)
        println(gL)

        print("\tCalculating ΦL: ")
        ΦL = Φ(currA, gL, L)
        println(ΦL)
        
        push!(gLs, gL)
        push!(ΦLs, ΦL)

        push!(As, nextA(ΦL, gL))
        println("\tNew A: $(As[end])")

        cond_val = abs(ΦL/Vₛ)
        println("\tabs(ΦL/Vₛ) = $cond_val")
        cond_val ≥ threshold || break
    end

    return As, ΦLs, gLs, x -> g(As[end], x), (g_x, x) -> Φ(As[end], g_x, x) 
end

result = simulate(1.97e-11, 2400, 8e6, 0.78, 1.5e-3, 0.5)
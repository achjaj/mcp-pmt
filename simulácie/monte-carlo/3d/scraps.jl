function atan2(x, y)
    x > 0 && y >= 0 && return atan(y/x)
    x < 0 && return atan(y/x) + π
    x > 0 && y < 0 && return atan(y/x) + 2π
    x == 0 && y > 0 && return π/2
    x == 0 && y < 0 && return 3π/2 
end

cart2spher(vec) = [norm(vec), atan2(vec[1:2]...), acos(vec[3]/norm(vec))]

function ddde(E0, NB, NP)
    D = E0/NB
    bins = [i*D for i in 0:1:NB]
    N = zeros(NB)
    Ns = 0
    Pn_hist = zeros(11)

    for p in 1:NP
        sey = see(E0, 0)
        l = length(sey)
        Ns += l
        Pn_hist[l + 1] += 1

        for se in sey
            if se <= bins[1]
                N[1] += 1
                continue
            end

            for i in 2:1:NB
                if bins[i - 1] < se <= bins[i]
                    N[i] += 1
                    break
                end
            end
        end
    end

    # normalization
    δe, δr, δ̃ts = getδs(E0, 0)
    δ = δe + δr + δ̃ts*(1 - δe - δr)

    mul = δ/(D*Ns)

    N .* mul, Pn_hist ./ Ns, bins
end
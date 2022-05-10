using Distributions
using Dates
using HDF5
using Base.Threads
using Mmap
using PlotlyJS
using StatsBase
using JSON

#= Notes:
    - energies T must be normed to T(0), where T = f(θ)

    - electron: array: [T₀, Tc, ϕ, θ, l₀, l, idx]
=#

lck = ReentrantLock()

⧷(A, B) = setdiff(A, B)

function safe_push!(arr, val)
    lock(lck) do
        push!(arr, val)
    end
end

function safe_append!(arr, val)
    lock(lck) do
        append!(arr, val)
    end
end

function safe_setindex!(A, X, inds...)
    lock(lck) do
        setindex!(A, X, inds...)
    end
end

function safe_getindex(A, inds...)
    val = nothing
    lock(lck) do
        val = getindex(A, inds...)
    end

    return val
end

function safe_println(msg)
    lock(lck) do
        println("[$(Threads.threadid())] $msg")
    end
end

sef(T, θ, β=0.6; a=0.62) = (T * sqrt(cos(θ)))^β * exp(a * (1 - cos(θ)) + β * (1 - T * sqrt(cos(θ))))                           # secondary emission function
vc(T₀::Float64, ϕ::Float64, q::Float64, E::Float64, m::Float64, d::Float64) =
    sqrt(2 * T₀ / m) * cos(ϕ) + q * E / m * sqrt(m / (2 * T₀)) * d / sin(ϕ), sqrt(2 * T₀ / m) * sin(ϕ)            # velocity at collision
lc(T₀::Float64, ϕ::Float64, q::Float64, E::Float64, d::Float64) = d / tan(ϕ) + (q * E * d^2) / (4 * T₀ * sin(ϕ)^2)       # point of collision

function θ(T₀::Float64, ϕ::Float64, q::Float64, E::Float64, m::Float64, d::Float64)                                      # collision angle
    vcx, vcy = vc(T₀, ϕ, q, E, m, d)
    acos(vcy / sqrt(vcx^2 + vcy^2))
end

function Tc(T₀::Float64, ϕ::Float64, q::Float64, E::Float64, m::Float64, d::Float64)                                     # collision energy
    vcx, vcy = vc(T₀, ϕ, q, E, m, d)
    1 / 2 * m * (vcx^2 + vcy^2)
end

function create_electron(T₀::Float64, ϕ::Float64, l::Float64, idx::Int)
    [T₀, -1.0, ϕ, -1.0, l, l, idx]
end

create_electron(parentT₀::Float64, l::Float64, idx::Int, rd::Rayleigh, cd::Cosine) = create_electron(parentT₀ * rand(rd), rand(cd), l, idx)

function create_step(group::HDF5.Group, index::Int, n::Int)
    step = create_dataset(group, "$index", datatype(Float64), dataspace(n, 7))
    step[1, 1] = 0
    return step
end

function simulate(E::Float64, d::Float64, L::Float64, initT₀::Float64, initϕ::Float64, initl₀::Float64; q=1.602e-19, ray_scale=1, μ=deg2rad(20), σ=deg2rad(10), βsmall=0.6, m=9.109e-31, max_step=-1, steps_path="steps.hdf5")
    rd = Rayleigh(ray_scale)
    cd = Cosine(μ, σ)

    wf = h5open(steps_path, "w")
    steps = create_group(wf, "steps")

    step_len = 1
    e_idx = 0
    step = create_step(steps, 1, step_len)
    step[1, :] = create_electron(initT₀, initϕ, initl₀, e_idx)
    e_idx += 1

    output = create_dataset(wf, "output", datatype(Float64), dataspace(1_000_000, 7))
    output_index = 1

    step_counter = 0
    while (max_step > -1 && step_counter < max_step) || (max_step <= -1 && step_len > 0)
        step_counter += 1
        println("Step $step_counter / $max_step; Input count: $step_len")

        electron_counts = Int[]
        T₀s = Float64[]
        ls = Float64[]
        left = Int[]
        @threads for i in 1:step_len
            GC.safepoint()

            params = safe_getindex(step, i, 1:2:3) # T₀, ϕ
            lcol = lc(params..., q, E, d) + safe_getindex(step, i, 6)
            safe_setindex!(step, lcol, i, 6)

            has_left = lcol > L
            if has_left
                safe_append!(left, i)

                safe_setindex!(output, safe_getindex(step, i, 1:7), output_index, 1:7)
                output_index += 1
            else
                Tcol = Tc(params..., q, E, m, d)
                θcol = θ(params..., q, E, m, d)
                safe_setindex!(step, Tcol, i, 2)
                safe_setindex!(step, θcol, i, 4)
                mean = sef(Tcol / params[1], θcol)
                safe_append!(electron_counts, rand(Poisson(mean)))
                safe_append!(T₀s, params[1])
                safe_append!(ls, lcol)
            end
        end

        stayed = collect(1:step_len) ⧷ left
        sl = length(stayed)

        next_step_len = sl + sum(electron_counts)            # No. of all electron in new step = no. of electrons in previous step + total no. of created electrons - no. of electrons that left the channel
        if next_step_len < 1
            break
        end

        next_step = create_step(steps, step_counter + 1, next_step_len)     # Create new step
        next_step[1:sl, :] = read(step)[stayed, 1:7]                        # Add electrons from previous step to new step
        next_step[1:sl, 1] = read(step)[stayed, 1] .* rand(rd, sl)          # Set T₀ of old electrons in new step
        next_step[1:sl, 3] = rand(cd, sl)                                   # Set ϕ of old electrons in new step

        # add new electrons to new step
        for (i, n) in enumerate(electron_counts)
            @threads for j in 1:n
                safe_setindex!(next_step,
                    create_electron(T₀s[i], ls[i], e_idx, rd, cd),
                    sl + sum(safe_getindex(electron_counts, 1:(i-1))) + j, 1:7
                )

                e_idx += 1
            end
        end

        step = next_step
        step_len = next_step_len
    end

    close(wf)
end

# =========================================================================================================================================

queue_len = 1_000_000
queue_top_ptr = 1
queue_bottom_ptr = 1

queue_file = open("queue", "w+")
queue = mmap(queue_file, Matrix{Float64}, (queue_len, 7))

output_file = open("output-$(now())", "w")
write(output_file, "[\n")

function add_to_queue(electron::Vector{Float64})
    global queue_bottom_ptr

    queue[queue_bottom_ptr, :] = electron
    Mmap.sync!(queue)
    global queue_bottom_ptr += 1
end

function add_to_output(electron::Vector)
    write(output_file, "\t$electron,\n")
    #flush(output_file)
end

function get_top(count=Threads.nthreads())
    global queue_top_ptr
    global queue_bottom_ptr

    rend = queue_top_ptr + count - 1
    range = queue_top_ptr:(rend <= queue_bottom_ptr ? rend : (queue_bottom_ptr - 1))
    global queue_top_ptr = range[end] + 1

    [[[j] for j in queue[i, :]] for i in range]     # converts array to array of arrays, so batch is array of arrays of arrays (electron is array of arrays)  
end


function simulate2(E::Float64, d::Float64, L::Float64, initT₀::Float64, initϕ::Float64, initl₀::Float64; q=1.602e-19, ray_scale=1, μ=deg2rad(20), σ=deg2rad(10), βsmall=0.6, m=9.109e-31, max_step=-1)
    rd = Rayleigh(ray_scale)
    cd = Cosine(μ, σ)

    e_idx = 0
    init_e = create_electron(initT₀, initϕ, initl₀, e_idx)
    e_idx = 1

    add_to_queue(init_e)

    while (queue_top_ptr == 1 && queue_top_ptr == queue_bottom_ptr) || (queue_top_ptr < queue_bottom_ptr)
        batch = get_top()
        left = Int[]

        println("Processing electrons: $(getindex.(batch, 7)) ($(queue_top_ptr - 1) / $(queue_bottom_ptr - 1) / $queue_len)")
        while length(left) != length(batch)
            @threads for (bi, electron) in collect(enumerate(batch))
                params = electron[1][end], electron[3][end]
                lcol = lc(params..., q, E, d) + electron[6][end]                            # calculate the point of collision
                append!(electron[6], lcol)                                                  # add the point of collision to electron data

                if bi in left
                    continue
                end

                if lcol > L                                                                 # test if electron leaves the channel
                    safe_println("Electron $(electron[end]) left the channel")
                    add_to_output(electron)                                                 # add the electron to output

                    safe_append!(left, bi)
                else
                    Tcol = Tc(params..., q, E, m, d)                                        # calculate the collision energy
                    θcol = θ(params..., q, E, m, d)                                         # calculate the collision angle

                    append!(electron[2], Tcol)
                    append!(electron[4], θcol)

                    mean = sef(Tcol / params[1], θcol)                                      # calculate the mean No. of new electrons
                    N = rand(Poisson(mean))                                                 # sample Poisson dist. for new No. of electrons

                    for i in 1:N
                        lock(lck) do                                                        # we need this lock because, we are accessing an modifying the e_idx
                            #println("[$(threadid())] Adding electron $e_idx from $(electron[end][end])")
                            add_to_queue(create_electron(params[1], lcol, e_idx, rd, cd))   # add new electron to queue
                            e_idx += 1
                        end
                    end

                    append!(electron[1], rand(rd) * params[1])                              # add new T₀ to the electron data
                    append!(electron[3], rand(cd))                                          # add new ϕ to the electron data
                end
            end
        end
    end
end

leaveTime(T₀::Float64, ϕ::Float64, q::Float64, E::Float64, m::Float64, L::Float64, l₀::Float64) = (sqrt(2 * T₀ * m * cos(ϕ)^2 + 2 * q * E * m * (L - l₀)) - sqrt(2 * T₀ * m) * cos(ϕ)) / (q * E)
v²(t::Float64, T₀::Float64, ϕ::Float64, m::Float64, q::Float64, E::Float64) = 2T₀ / m + q^2 * E^2 / m^2 * t^2 + sqrt(2 * T₀ / m) * cos(ϕ) * q * E / m * t
Tleave(T₀::Float64, ϕ::Float64, m::Float64, q::Float64, E::Float64, L::Float64, l₀::Float64) = 1 / 2 * m * v²(leaveTime(T₀, ϕ, q, E, m, L, l₀), T₀, ϕ, m, q, E)
leaveAng(T₀::Float64, ϕ::Float64, m::Float64, q::Float64, E::Float64, L::Float64, l₀::Float64) = acos((sqrt(2T₀ / m) * sin(ϕ)) / sqrt(v²(leaveTime(T₀, ϕ, q, E, m, L, l₀), T₀, ϕ, m, q, E)))

d = 10e-6  # m
L = 40 * d   # m
E = 4000 / L # V/m
eV = 1.602e-19 # J

simulate2(E, d, 20 * d, 100.0, deg2rad(45), 0.0, m=0.511)

skip(output_file, -2)
write(output_file, "\n]")
close(output_file)
close(queue_file)

# ==============================================================================================
function analyse(path)
    file = open(path, "r")
    electrons = JSON.parse(line)
    file.close()

    Tls = [Tleave(electron[1][end], electron[3][end], 0.511, 1.602e-19, 4000 / (40 * 10e-6), 40 * 10e-6, electron[6][end-1]) for electron in electrons]
    TlsHist = fit(Histogram, Tls, 0:1:10e3)
    TlsHP = scatter(x=TlsHist.edges[1], y=TlsHist.weights)

    T₀s = vcat(getindex.(electrons, 1)...)
    T₀sHist = fit(Histogram, T₀s, 0:100:10e3)
    T₀sHP = scatter(x=T₀sHist.edges[1], y=T₀sHist.weights)

    Tcs = filter(x -> x != -1, vcat(getindex.(electrons, 2)...))
    TcsHist = fit(Histogram, Tcs, 0:100:10e3)
    TcsHP = scatter(x=TcsHist.edges[1], y=TcsHist.weights)

    ϕs = vcat(getindex.(electrons, 3)...)
    ϕsHist = fit(Histogram, ϕs, 0:1e-3:maximum(ϕs))
    ϕsHP = scatter(x=ϕsHist.edges[1], y=ϕsHist.weights)

    θs = vcat(getindex.(electrons, 4)...)
    θsHist = fit(Histogram, θs, 1:1e-3:1.4)
    θsHP = scatter(x=θsHist.edges[1], y=θsHist.weights)

    las = rad2deg.([leaveAng(electron[1][end], electron[3][end], 0.511, 1.602e-19, 4000 / (40 * 10e-6), 40 * 10e-6, electron[6][end-1]) for electron in electrons])
    lasHist = fit(Histogram, las, rad2deg(1):rad2deg(1e-3):maximum(las))
    lasHP = scatter(x=lasHist.edges[1], y=lasHist.weights)

    dbc = vcat([diff(electron[end-1]) for electron in electrons]...)
    dbcHist = fit(Histogram, dbc, 1e-5:1e-7:maximum(dbc))
    dbcHP = scatter(x=dbcHist.edges[1], y=dbcHist.weights)
end

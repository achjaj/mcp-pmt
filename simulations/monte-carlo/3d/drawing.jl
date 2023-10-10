using GLMakie
using GeometryBasics: Cylinder, Point3f0
using Printf

rotMat(β) = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]

# creates plot with a cylinder, i. e. our channel
function preparePlot()
    #creates empty plot
    plt = Figure()
    ax = Axis3(plt[1, 1], aspect = (1, 3, 1), title = "V = $(@sprintf "%.0e" V) V" #=limits = (nothing, nothing, 0, L, nothing, nothing)=#)

    # draw the cylinder in red color and transparency
    mesh!(ax, Cylinder{3, Float32}(Point3f0(0.0, 0.0, R), Point3f0(0.0, L, R), R), color = (:red, 0.3), transparency=true)
    
    return ax, plt
end


function drawSteps(steps, res = 100)
    colors = [:blue, :red]
    ax, plt = preparePlot()

    r0 = zeros(3)
    for (i, s) in enumerate(steps)
        last = i == length(steps)
        d̄ = r0[[1, 3]]
        d = norm(d̄)

        if d > 2R
            @warn "Step $i: d > 2R ($d̄, $d)"
            d = 2R
        end

        ratio = d/R/2
        abs(ratio) > 1 && @error "Step $i: ratio > 1 ($ratio)"

        sgn = r0[1] > 0 ? -1 : 1
        β = 2*asin(ratio)*sgn
        
        mat = rotMat(β)

        times = range(0, s[2], length = res)
        points = [mat*r(s[1], t)+r0 for t in times]

        last && filter!(point -> point[2] <= L, points)

        c = last ? :lightgreen : colors[(i % 2) + 1]
        lines!(ax, hcat(points...), color = c)

        r0 = points[end]
    end

    plt, ax
end

function drawPositions(positions)
    ax, plt = preparePlot()
    lines!(ax, hcat(positions...))

    plt, ax
end

function plotT₀(steps)
    plt = Figure()
    ax = Axis(plt[1, 1], xlabel = "Step", ylabel = "T₀ [eV]", subtitle = "V = $(@sprintf "%.0e" V) V", title = "Energy at the start of each step")

    v̄s = getindex.(steps, 1)
    vs = norm.(v̄s)
    Ts = [0.5*m*v^2 /q for v in vs]

    scatterlines!(ax, Ts)

    plt
end

function sVstc(steps)
    plt = Figure()
    ax = Axis(plt[1, 1], xlabel = "t [s]", ylabel = "y [m]", subtitle = "V = $(@sprintf "%.0e" V) V", title = "Position along the channel vs. time")

    tcs = [0; getindex.(steps, 2)]
    ss = [0; getindex.(steps, 3)]
    scatterlines!(ax, cumsum(tcs), ss)

    plt
end
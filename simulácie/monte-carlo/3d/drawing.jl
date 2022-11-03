using GLMakie
using GeometryBasics: Cylinder, Point3f0

rotMat(β) = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]

# creates plot with a cylinder, i. e. our channel
function preparePlot()
    #creates empty plot
    plt = Figure()
    ax = Axis3(plt[1, 1], aspect = (1, 3, 1), #=limits = (nothing, nothing, 0, L, nothing, nothing)=#)

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

        ratio = norm(d)/R/2
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
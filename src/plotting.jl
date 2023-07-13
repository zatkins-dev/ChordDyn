export LatexTicks, dt_color, trajectory!, trajectory3!, arrow3d!, plot, plot_pendulum_on_tonnetz, animate_pendulum, draw_pendulum!
import CairoMakie.plot

## Plotting Utilities
struct LatexTicks
    ticks
end

LatexTicks() = LatexTicks(Makie.MakieCore.automatic)

function Makie.get_ticks(t::LatexTicks, scale, formatter, vmin, vmax)
    values, labels = Makie.get_ticks(t.ticks, scale, formatter, vmin, vmax)
    values, map(latexstring, labels)
end

function dt_color(t; colormap=:seaborn_deep)
    dts = round.(t[2:end] - t[1:end-1]; sigdigits=4)
    bins = unique(sort(dts))
    push!(dts, dts[end])
    bbins = [0, bins..., 1]
    bounds = [(bbins[i] / 2 + bbins[i-1] / 2 for i in 2:length(bbins)-1)..., bbins[end-1] + bbins[end-2] / 2]
    indices = [findfirst(x -> x == dt, bins) for dt in dts]
    values = [(0:length(bins))...] / length(bins)
    middles = round.(((1:length(bins)) .- 0.5) / length(bins); digits=3)
    colors = map(i -> middles[i], indices)
    labels = ["$(round(v, sigdigits=2))" for v in bins]
    cm = cgrad(colormap, values; categorical=true)
    return cm, middles, labels, colors
end

function trajectory!(X, V; colors=1:size(X, 2), colormap=:batlow, label=undef, arrowsize=32, linewidth=4)
    arrows!(X[1, :], X[2, :], V[1, :], V[2, :], arrowsize=arrowsize, arrowcolor=colors, linecolor=colors, linewidth=linewidth, colormap=colormap, label=label)
end

function arrow3d!(x, y, z, u, v, w; as=0.1, lc=:black, lw=0.4)
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x, y, z, u, v, w, i) in zip(x, y, z, u, v, w, 1:length(x))
        nv = sqrt(u^2 + v^2 + w^2)
        v1, v2 = -[u, v, w] / nv, nullspace(adjoint([u, v, w]))[:, 1]
        v4 = (3 * v1 + v2) / 3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2 * (v4' * v2) * v2
        (as < 0) && (nv = nv0)
        v4, v5 = -as * nv * v4, -as * nv * v5
        color = lc
        if lc isa AbstractArray{<:Any,1} || lc isa AbstractRange{<:Number}
            color = lc[i]
        end
        lines!([x, x + u], [y, y + v], [z, z + w]; transparency=true, color=color, linewidth=lw, label=false)
        lines!([x + u, x + u - v5[1]], [y + v, y + v - v5[2]], [z + w, z + w - v5[3]]; transparency=true, color=color, linewidth=lw, label=false)
        lines!([x + u, x + u - v4[1]], [y + v, y + v - v4[2]], [z + w, z + w - v4[3]]; transparency=true, color=color, linewidth=lw, label=false)
    end
end

function trajectory3!(X, V; colors=:black)
    arrow3d!(X[1, :], X[2, :], X[3, :], V[1, :], V[2, :], V[3, :]; lc=colors)
end

## Plot functions

function plot(t::TriangularLattice; xres=1920)
    aspect = xlength(t) / ylength(t)
    f = Figure(resolution=(xres, xres / aspect))
    ax = Axis(f[1, 1]; aspect=aspect / 1)
    scatter!(points(t); markersize=12, label=nothing)
    xs = xlims(t)
    ys = ylims(t)
    resize_to_layout!(f)
    return f, ax
end

eharmonic(n::Int) = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "χ", "ε"][mod(n, 12)+1]
midi(n::Int) = "$n"
octave(n::Int) = "$(n ÷ 12)"
note(n::Int) = ["C", "C♯", "D", "D♯", "E", "F", "F♯", "G", "G♯", "A", "A♯", 'B'][mod(n, 12)+1]

vertexlabel_fns = (:eharmonic, :midi, :octave, :note)

function get_triangles(t::TriangularLattice, color=:blue)
    xs = xlims(t)
    ys = ylims(t)
    if color == :blue
        centers = points(t) .+ Vec2f(tri_x(t) / 2, 0)
    else
        centers = points(t) .- Vec2f(tri_x(t) / 2, 0)
    end
    tri_points = [Point2f.(p) for p in closest_points.(Ref(t), centers)]
    triangles = Polygon[]
    split_triangles = Polygon[]
    for pts in tri_points
        at_top = Vector{Bool}(undef, 3)
        at_bottom = Vector{Bool}(undef, 3)
        at_left = Vector{Bool}(undef, 3)
        at_right = Vector{Bool}(undef, 3)
        for (i, p) in enumerate(pts)
            at_top[i] = abs(p[2] + tri_y(t) / 2 - ys[2]) < 1e-6
            at_bottom[i] = abs(p[2] - ys[1]) < 1e-6
            at_left[i] = abs(p[1] - xs[1]) < 1e-6
            at_right[i] = abs(p[1] + tri_x(t) - xs[2]) < 1e-6
        end

        if any(at_right)
            for i in eachindex(pts)
                if at_left[i]
                    pts[i] += Vec2f(xs[2], 0)
                end
            end
        end

        if !(any(at_top) && any(at_bottom))
            push!(triangles, Polygon(pts))
        elseif any(at_top) && any(at_bottom)
            last_p = pts[.!at_top.&&.!at_bottom][1]
            if last_p[2] > ys[2] / 2
                for i in eachindex(pts)
                    if at_bottom[i]
                        pts[i] = Vec2f(pts[i][1], ys[2])
                    end
                end
                push!(triangles, Polygon(pts))
            else
                top_tri = deepcopy(pts)
                for i in eachindex(top_tri)
                    if !at_top[i]
                        top_tri[i] = Point2f(top_tri[i][1], ys[2])
                    end
                end
                push!(split_triangles, Polygon([top_tri[at_bottom][1], top_tri[at_top][1], top_tri[.!at_bottom.&&.!at_top][1]]))
                bottom_tri = deepcopy(pts)
                for i in eachindex(bottom_tri)
                    if at_top[i]
                        bottom_tri[i] = Point2f(bottom_tri[i][1], ys[1])
                    end
                end
                push!(split_triangles, Polygon([bottom_tri[at_bottom][1], bottom_tri[.!at_bottom.&&.!at_top][1], bottom_tri[at_top][1]]))
            end
        elseif any(at_top) && !any(at_bottom) && !(any(at_left) && any(at_right))
            push!(triangles, Polygon(pts))
        end
    end
    return triangles, split_triangles
end

function CairoMakie.plot(tonnetz::Tonnetz; xres=1920, vertexlabels=:eharmonic)
    @assert vertexlabels in vertexlabel_fns "vertexlabels must be one of $(vertexlabel_fns)"

    aspect = xlength(tonnetz) / ylength(tonnetz)
    f = Figure(resolution=(xres, xres / aspect), backgroundcolor=:transparent)
    ax = Axis(f[1, 1], xgridvisible=false, ygridvisible=false, xticksvisible=false, yticksvisible=false, backgroundcolor=:transparent)
    hidespines!(ax)
    hidedecorations!(ax)

    for color in [:blue, :red]
        triangles, split_triangles = get_triangles(lattice(tonnetz), color)
        for t in triangles
            poly!(t, color=(color, 0.2), label=nothing, strokewidth=0, transparent=true)
            c = coordinates(t)
            lines!([c..., c[1]]; color=:black, linewidth=2, label=nothing, transparent=true)
        end
        for t in split_triangles
            poly!(t; color=(color, 0.2), label=nothing, strokewidth=0, transparent=true)
            lines!(coordinates(t); color=:black, linewidth=2, label=nothing, transparent=true)
        end
    end
    split_triangles = Polygon[]

    markers = eval(:($vertexlabels.(notes($tonnetz))))
    longest_marker = maximum(length.(markers))
    textsize = 32 - 4 * (longest_marker - 1)
    scatter!(points(tonnetz); marker=:circle, color=:black, markersize=64, label=nothing, transparent=false)
    scatter!(points(tonnetz); marker=:circle, color=:white, markersize=56, label=nothing, transparent=false)
    text!(points(tonnetz); text="" .* (markers), font="ComicCodeLigatures NF", align=(:center, :center), color=:black, textsize=textsize, transparent=false)
    left_col = points(tonnetz)[1:12] .+ Vec2f(xlength(tonnetz), 0)
    scatter!(left_col; marker=:circle, color=:black, markersize=64, label=nothing, transparent=false)
    scatter!(left_col; marker=:circle, color=:white, markersize=56, label=nothing, transparent=false)
    text!(left_col; text="" .* (markers[1:12]), font="ComicCodeLigatures NF", align=(:center, :center), color=:black, textsize=textsize, transparent=false)
    top_row = points(tonnetz)[1:24:end] .+ Vec2f(0, ylength(tonnetz))
    scatter!(top_row; marker=:circle, color=:black, markersize=64, label=nothing, transparent=false)
    scatter!(top_row; marker=:circle, color=:white, markersize=56, label=nothing, transparent=false)
    text!(top_row; text="" .* (markers[1:24:end]), font="ComicCodeLigatures NF", align=(:center, :center), color=:black, textsize=textsize, transparent=false)
    tl = points(tonnetz)[1] .+ Vec2f(xlength(tonnetz), ylength(tonnetz))
    scatter!(tl; marker=:circle, color=:black, markersize=64, label=nothing, transparent=false)
    scatter!(tl; marker=:circle, color=:white, markersize=56, label=nothing, transparent=false)
    text!(tl; text="" * markers[1], font="ComicCodeLigatures NF", align=(:center, :center), color=:black, textsize=textsize, transparent=false)
    resize_to_layout!(f)
    return f, ax
end

function plot_pendulum_on_tonnetz(tonnetz, X, t; arrowsize=24, linewidth=4)
    fig, ax = plot(tonnetz)
    pts = project_onto.(Ref(tonnetz.t), X[1, :], X[2, :])
    θ₁ = [p[1] for p in pts]
    θ₂ = [p[2] for p in pts]
    xs = hcat(θ₁, θ₂)'
    vx = θ₁[2:end] - θ₁[1:end-1]
    vy = θ₂[2:end] - θ₂[1:end-1]
    for i in eachindex(vx)
        if abs(vx[i]) > 0.75 * xlength(tonnetz)
            vx[i] = -sign(vx[i]) * (xlength(tonnetz) - abs(vx[i]))
        end
        if abs(vy[i]) > 0.75 * ylength(tonnetz)
            vy[i] = -sign(vy[i]) * (ylength(tonnetz) - abs(vy[i]))
        end
    end
    V = hcat(vx, vy)'
    trajectory!(xs[:, 1:end-1], V; colors=:black, arrowsize=arrowsize + 3, linewidth=linewidth + 1.5)
    trajectory!(xs[:, 1:end-1], V, colors=t[1:end-1]; arrowsize=arrowsize, linewidth=linewidth)
    return fig, ax
end

function draw_pendulum!(dp::DoublePendulum, θ₁, θ₂, t)
    x1, y1 = dp.ℓ₁ * sin(θ₁), -dp.ℓ₁ * cos(θ₁)
    x2, y2 = x1 + dp.ℓ₂ * sin(θ₂), y1 - dp.ℓ₂ * cos(θ₂)

    lines!([0, x1], [0, y1], color=:black, linewidth=10)
    lines!([0, x1], [0, y1]; color=:red, linewidth=6)
    lines!([x1, x2], [y1, y2], color=:black, linewidth=10)
    lines!([x1, x2], [y1, y2]; color=:red, linewidth=6)
    scatter!(0, 0, color=:black, markersize=20)

    scatter!(x1, y1; color=:black, markersize=96)
    scatter!(x1, y1; color=:white, markersize=88)
    text!(x1, y1; text=L"m_1", align=(:center, :center), color=:black, textsize=40, markerspace=:pixel)
    scatter!(x2, y2; color=:black, markersize=96 * dp.m₂ / dp.m₁)
    scatter!(x2, y2; color=:white, markersize=88 * dp.m₂ / dp.m₁)
    text!(x2, y2; text=L"m_2", align=(:center, :center), color=:black, textsize=40, markerspace=:pixel)
    text!(0.5, 1, align=(:center, :top), text=latexstring("t = $(round(t, digits = 1))"), textsize=48, space=:relative)
end

function animate_pendulum(dp::DoublePendulum, θ₁, θ₂, t; output="output/pendulum_anim.mp4", framerate=nothing)
    tstep = Observable(1)
    dt = t[2] - t[1]

    x1, y1 = @lift(dp.ℓ₁ * sin(θ₁[$tstep])), @lift(-dp.ℓ₁ * cos(θ₁[$tstep]))
    x2, y2 = @lift($x1 + dp.ℓ₂ * sin(θ₂[$tstep])), @lift($y1 - dp.ℓ₂ * cos(θ₂[$tstep]))
    f = Figure(resolution=(1200, 1200), backgroundcolor=:white)
    ax = Axis(f[1, 1],
        aspect=1,
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        backgroundcolor=:white,
        limits=(-dp.ℓ₁ - dp.ℓ₂ - 0.1, dp.ℓ₁ + dp.ℓ₂ + 0.1, -dp.ℓ₁ - dp.ℓ₂ - 0.1, dp.ℓ₁ + dp.ℓ₂ + 0.1),
    )
    hidespines!(ax)
    hidedecorations!(ax)

    lines!(@lift([0, $x1]), @lift([0, $y1]), color=:black, linewidth=10)
    lines!(@lift([0, $x1]), @lift([0, $y1]); color=:red, linewidth=6)
    lines!(@lift([$x1, $x2]), @lift([$y1, $y2]), color=:black, linewidth=10)
    lines!(@lift([$x1, $x2]), @lift([$y1, $y2]); color=:red, linewidth=6)
    scatter!(0, 0, color=:black, markersize=20)

    scatter!(x1, y1; color=:black, markersize=96)
    scatter!(x1, y1; color=:white, markersize=88)
    text!(x1, y1; text=L"m_1", align=(:center, :center), color=:black, textsize=40, markerspace=:pixel)
    scatter!(x2, y2; color=:black, markersize=96 * dp.m₂ / dp.m₁)
    scatter!(x2, y2; color=:white, markersize=88 * dp.m₂ / dp.m₁)
    text!(x2, y2; text=L"m_2", align=(:center, :center), color=:black, textsize=40, markerspace=:pixel)
    text!(0.5, 1, align=(:center, :top), text=@lift(latexstring("t = $(round(t[$tstep], digits = 1))")), textsize=48, space=:relative)

    itr = nothing
    if isnothing(framerate)
        framerate = round(Int, 1 / dt)
        itr = eachindex(t)
    else
        itr = round.(range(1, length(t), length=round(Int, (t[end] - t[1]) * framerate)))
    end
    record(f, output, itr; framerate=framerate) do frame
        tstep[] = frame
    end
end
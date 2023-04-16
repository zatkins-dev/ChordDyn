export LatexTicks, dt_color, trajectory!, trajectory3!, arrow3d!, plot, plot_pendulum_on_tonnetz, animate_pendulum, note_marker

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

function trajectory!(X, V; colors=1:size(X, 2), colormap=:seaborn_deep, label=undef)
    arrows!(X[1, :], X[2, :], V[1, :], V[2, :], arrowsize=8, arrowcolor=colors, linecolor=colors, linewidth=1, colormap=colormap, label=label)
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

note_marker(n::Int) = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'χ', 'ε'][mod(n, 12)+1]

function plot(tonnetz::Tonnetz; xres=1920)
    markers = note_marker.(notes(tonnetz))
    aspect = xlength(tonnetz) / ylength(tonnetz)
    f = Figure(resolution=(xres, xres / aspect))
    ax = Axis(f[1, 1]; aspect=aspect)
    p = scatter!(points(tonnetz); marker=markers, markersize=32, label=nothing)
    xs = xlims(tonnetz)
    ys = ylims(tonnetz)
    domain = Rect2f((0, 0), (xlength(tonnetz), ylength(tonnetz)))
    recmesh = GeometryBasics.mesh(domain)
    Makie.mesh!(recmesh; color=(:orange, 0.1))
    xlims!(xs[1] - 0.1, xs[2] + 0.1)
    ylims!(ys[1] - 0.1, ys[2] + 0.1)
    resize_to_layout!(f)
    return f, ax
end

function plot_pendulum_on_tonnetz(tonnetz, X, t)
    fig, ax = plot(tonnetz)
    pts = project_onto.(Ref(tonnetz.t), X[1, :], X[2, :])
    θ₁ = [p[1] for p in pts]
    θ₂ = [p[2] for p in pts]
    xs = hcat(θ₁, θ₂)'
    V = xs[:, 2:end] - xs[:, 1:end-1]
    trajectory!(xs[:, 1:end-1], V, colors=t[1:end-1])
    save("output/pendulum_on_tonnetz.png", fig, px_per_unit=2)
    display(fig)
end

function animate_pendulum(dp::DoublePendulum, θ₁, θ₂, t; output="output/pendulum_anim.mp4", framerate=nothing)
    tstep = Observable(1)
    dt = t[2] - t[1]

    x1, y1 = @lift(dp.ℓ₁ * sin(θ₁[$tstep])), @lift(-dp.ℓ₁ * cos(θ₁[$tstep]))
    x2, y2 = @lift($x1 + dp.ℓ₂ * sin(θ₂[$tstep])), @lift($y1 - dp.ℓ₂ * cos(θ₂[$tstep]))

    fig, ax = lines(@lift([0, $x1]), @lift([0, $y1]), color=:blue, linewidth=4,
        axis=(title=@lift("t = $(round(t[$tstep], digits = 1))"),))

    lines!(@lift([$x1, $x2]), @lift([$y1, $y2]); color=:red, linewidth=4)
    scatter!(x1, y1; color=:blue, markersize=50)
    scatter!(x2, y2; color=:blue, markersize=50 * dp.m₂ / dp.m₁)

    total_length = dp.ℓ₁ + dp.ℓ₂
    pad = 0.1 * total_length
    limits!(ax, -total_length - pad, total_length + pad, -total_length - pad, total_length + pad)

    itr = nothing
    if isnothing(framerate)
        framerate = round(Int, 1 / dt)
        itr = eachindex(t)
    else
        itr = round.(range(1, length(t), length=round(Int, (t[end] - t[1]) * framerate)))
    end
    record(fig, output, itr;
        framerate=framerate) do frame
        tstep[] = frame
    end
end
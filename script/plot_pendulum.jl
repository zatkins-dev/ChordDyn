import ChordDyn, Ordinary
import CairoMakie
import LaTeXStrings


let
    dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]
    dt = 1e-3

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=100, Δt=dt, returnvectors=false)

    transient = findfirst(t .>= 90)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]

    θ₁ = mod2pi.(X[1, :])
    θ₂ = mod2pi.(X[2, :])

    x1, y1 = dp.ℓ₁ * sin(θ₁[1]), -dp.ℓ₁ * cos(θ₁[1])
    x2, y2 = x1 + dp.ℓ₂ * sin(θ₂[1]), y1 - dp.ℓ₂ * cos(θ₂[1])
    x1_n, y1_n = dp.ℓ₁ * sin(θ₁[1] + X[3, 1]), -dp.ℓ₁ * cos(θ₁[1] + X[3, 1])
    x2_n, y2_n = x1_n + dp.ℓ₂ * sin(θ₂[1] + X[4, 1]), y1_n - dp.ℓ₂ * cos(θ₂[1] + X[4, 1])

    fig = Figure(resolution=(1080, 1080), backgroundcolor=:transparent)
    ax = Axis(fig[1, 1], aspect=1, xlabel=L"x", ylabel=L"y", backgroundcolor=:transparent,
        limits=(-dp.ℓ₁ - dp.ℓ₂ - 0.1, dp.ℓ₁ + dp.ℓ₂ + 0.1, -dp.ℓ₁ - dp.ℓ₂ - 0.1, dp.ℓ₁ + dp.ℓ₂ + 0.1),
    )
    arrows!([x1], [y1], [x1_n - x1], [y1_n - y1]; arrowsize=25, arrowcolor=:red, linestyle=:dash, linecolor=:red, linewidth=6)
    arrows!([x2], [y2], [x2_n - x2], [y2_n - y2]; arrowsize=25, arrowcolor=:red, linestyle=:dash, linecolor=:red, linewidth=6)
    draw_pendulum!(dp, θ₁[1], θ₂[1], t[1])
    save("output/pendulum.png", fig, px_per_unit=2)
end

let
    dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]
    dt = 1e-3

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=300, Δt=dt, returnvectors=false)

    transient = findfirst(t .>= 60)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]
    xlim = extrema(X[1, :])
    ylim = extrema(X[2, :])
    aspect = (xlim[2] - xlim[1]) / (ylim[2] - ylim[1])

    fig = Figure(resolution=(800 * aspect + 50, 800), fontsize=20, backgroundcolor=:transparent)
    ax = Axis(fig[1, 1], aspect=aspect, xlabel=L"\theta_1", ylabel=L"\theta_2", xticks=LatexTicks(MultiplesTicks(7, pi, "π")), yticks=LatexTicks(MultiplesTicks(7, pi, "π")), backgroundcolor=:transparent)
    V = X[1:2, 2:end] - X[1:2, 1:end-1]
    # for i in eachindex(V[1, :])
    #     if abs(V[1, i]) > 0.75 * 2π
    #         V[1, i] = -sign(V[1, i]) * (2π - abs(V[1, i]))
    #     end
    #     if abs(V[2, i]) > 0.75 * 2π
    #         V[2, i] = -sign(V[2, i]) * (2π - abs(V[2, i]))
    #     end
    # end
    p = trajectory!(X[1:2, 1:end-1], V; linewidth=1, arrowsize=2)
    Colorbar(fig[1, 2], label=L"\text{Time}", colormap=:batlow, limits=extrema(t), ticks=LatexTicks())
    save("output/pendulum_trajectory.png", fig, px_per_unit=2)
end

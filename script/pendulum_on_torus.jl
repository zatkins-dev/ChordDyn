using Ordinary
using ChordDyn
using GLMakie
using ColorSchemes

let
    dp = DoublePendulum(m₁=0.1, m₂=0.2, ℓ₁=0.5, ℓ₂=1)
    x₀ = [π, π - 0.3, 0, 0]
    dt = AdaptTS(errₘₐₓ=1e-3)
    X, V, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=10, Δt=dt, returnvectors=true)
    X[1, :] = mod2pi.(X[1, :])
    X[2, :] = mod2pi.(X[2, :])

    mod12(θ) = mod2pi(θ) / π * 6
    mod12to2pi(n) = mod(n, 12) / 6 * π

    torus(θ₁, θ₂; a=5, c=10) = [
        (c + a * cos(θ₂)) * cos(θ₁),
        (c + a * cos(θ₂)) * sin(θ₁),
        a * sin(θ₂)
    ]

    a = 5
    c = 10

    # plot torus
    Θ = ϕ = range(-π, π, 300)
    x = [torus(u, v; a=a, c=c)[1] for u in Θ, v in ϕ]
    y = [torus(u, v; a=a, c=c)[2] for u in Θ, v in ϕ]
    z = [torus(u, v; a=a, c=c)[3] for u in Θ, v in ϕ]
    fig = surface(x, y, z, colormap=[(:grey)],
        lightposition=Vec3f(0, 0, 0.8), ambient=Vec3f(0.6, 0.6, 0.6),
        backlight=2.0f0)

    # plot trajectory
    tor = reduce(hcat, torus.(X[1, :], X[2, :]; a=a + 0.25, c=c))
    tor_V = tor[:, 2:end] - tor[:, 1:end-1]
    colors = get.(Ref(ColorSchemes.jet), range(0, 1, length=length(tor_V[1, 1:end])))
    trajectory3!(tor[:, 1:end-1], tor_V, colors=colors)
    # scatter!(tor[1, :], tor[2, :], tor[3, :]; color=1:length(tor[1, :]), colormap=:jet)

    # plot tonnetz
    tonnetz = Tonnetz()
    apx, apy = angular_points(tonnetz)
    tor_tonnetz = reduce(hcat, torus.(apx, apy; a=a + 0.25, c=c))
    scatter!(tor_tonnetz[1, :], tor_tonnetz[2, :], tor_tonnetz[3, :]; marker=note_marker.(notes(tonnetz)), markersize=32, label=nothing)

    display(fig)
end
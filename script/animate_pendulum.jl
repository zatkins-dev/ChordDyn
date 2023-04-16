using Ordinary, ChordDyn

let
    dp = DoublePendulum(m₁=2.0, m₂=1.0, ℓ₁=2.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=70, Δt=1e-3)

    transient = findfirst(t .>= 60)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]

    θ₁ = mod2pi.(X[1, :])
    θ₂ = mod2pi.(X[2, :])

    animate_pendulum(dp, θ₁, θ₂, t; framerate=120)
end
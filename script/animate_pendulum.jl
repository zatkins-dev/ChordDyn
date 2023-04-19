using Ordinary, ChordDyn

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

    animate_pendulum(dp, θ₁, θ₂, t; framerate=120)
end
using Ordinary, ChordDyn

let
  dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
  x₀ = [π, π - 0.3, 0, 0]
  dt = 1e-3

  X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=180, Δt=dt, returnvectors=false)

  transient = findfirst(t .>= 60)#Int(5 / dt)
  X = X[:, transient+1:end]
  t = t[transient+1:end]


  animate_tonnetz(dp, X, t; framerate=120, vertexlabels=:note)
end

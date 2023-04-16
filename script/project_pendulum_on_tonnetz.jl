using Ordinary, ChordDyn
using DelimitedFiles

let
    dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]
    dt = 1e-3

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=300, Δt=dt, returnvectors=false)

    transient = findfirst(t .>= 60)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]
    tonnetz = Tonnetz()
    X[1:2, :] = X[2:-1:1, :]
    # Xhat, that = X, t
    # SECTIONING WITH AROUND 0.25 WORKS QUITE WELL
    s = TemporalSection(0.125)
    Xhat, that = section(s, X, t; interpolate=true)
    plot_pendulum_on_tonnetz(tonnetz, Xhat, that)
    chords, durations = trajectory_to_chord_progression(tonnetz, Xhat[2, :], Xhat[1, :], that)

    writedlm("output/chords.txt", chords, " ")
    writedlm("output/durations.txt", 2 * durations, " ")
end
using Ordinary, ChordDyn
using DelimitedFiles

let
    dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]
    dt = 1e-3

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=100, Δt=dt, returnvectors=false)

    transient = findfirst(t .>= 90)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]
    tonnetz = Tonnetz()
    chords, durations = trajectory_to_chord_progression(tonnetz, X[1, :], X[2, :], t)

    writedlm("output/chords.txt", chords, " ")
    writedlm("output/durations.txt", 2 * durations, " ")


    # SECTIONING WITH AROUND 0.125 WORKS QUITE WELL
    s = TemporalSection(0.125)
    Xhat, that = section(s, X, t; interpolate=false)
    plot_pendulum_on_tonnetz(tonnetz, Xhat, that)
    @show Xhat[1:2, end-4:end]'
    @show that[end-4:end]
    section_chords, section_durations = trajectory_to_chord_progression(tonnetz, Xhat[1, :], Xhat[2, :], that)

    writedlm("output/chords-section.txt", section_chords, " ")
    writedlm("output/durations-section.txt", 2 * section_durations, " ")
end
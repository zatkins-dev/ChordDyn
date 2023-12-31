using Ordinary, ChordDyn
using DelimitedFiles
using CairoMakie

let
    dp = DoublePendulum(m₁=1.0, m₂=1.0, ℓ₁=1.0, ℓ₂=1.0)
    x₀ = [π, π - 0.3, 0, 0]
    dt = 1e-3

    X, t = solve(rk4, ∂ₜ(dp); x₀=x₀, t₀=0, tₑ=180, Δt=dt, returnvectors=false)

    transient = findfirst(t .>= 60)#Int(5 / dt)
    X = X[:, transient+1:end]
    t = t[transient+1:end]

    tonnetz = Tonnetz()
    s = TemporalSection(0.125)
    Xhat, that = section(s, X, t; interpolate=false)
    fig, ax = plot_pendulum_on_tonnetz(tonnetz, X, t, arrows=Xhat)
    save("output/pendulum_on_tonnetz_octave.png", fig, px_per_unit=2)
    chords, durations = trajectory_to_chord_progression(tonnetz, X[1, :], X[2, :], t)

    writedlm("output/chords.txt", chords, " ")
    writedlm("output/durations.txt", durations, " ")


    # # SECTIONING WITH AROUND 0.125 WORKS QUITE WELL
    # fig, ax = plot_pendulum_on_tonnetz(tonnetz, Xhat, that)
    # save("output/pendulum_on_tonnetz-section.png", fig, px_per_unit=2)
    # display(fig)
    # @show Xhat[1:2, end-4:end]'
    # @show that[end-4:end]
    section_chords, section_durations = trajectory_to_chord_progression(tonnetz, Xhat[1, :], Xhat[2, :], that)

    writedlm("output/chords-section.txt", section_chords, " ")
    writedlm("output/durations-section.txt", section_durations, " ")
end

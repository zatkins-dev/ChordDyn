export project_onto, trajectory_to_chord_progression

function project_onto(t::TriangularLattice, θ1, θ2)
    """Convert θ1, θ2 ∈ [0, 2π]² to t.xlims × t.ylims"""
    local x = rescale_to_lims(mod(θ1, 16π) / 16π, xlims(t))
    local y = rescale_to_lims(mod(θ2, 12π) / 12π, ylims(t))
    x, y
end

function closest_points(t::TriangularLattice, p)
    dims = Vec2f(xlength(t), ylength(t))
    p_wrapped = Vec2f(mod.(p, dims))
    distances = [periodic_distance(point, p_wrapped; dims=dims) for point in points(t)]
    closest = partialsortperm(norm.(distances), 1:3)
    return points(t)[closest]
end

function closest_points(tonnetz::Tonnetz, p)
    closest_points(lattice(tonnetz), p)
end

function project_onto(tonnetz::Tonnetz, θ1, θ2, notemode=false)
    """Convert θ1, θ2 ∈ [0, 2π]² to triplet of notes in [0,...,12]"""
    # use θ1 to compute x
    xy = Vec2f(project_onto(lattice(tonnetz), Float64(θ1), Float64(θ2)))
    dims = Vec2f(xlength(tonnetz), ylength(tonnetz))
    distances = [periodic_distance(p, xy; dims=dims) for p in points(tonnetz)]
    if notemode
        closest = argmin(norm.(distances))
        return notes(tonnetz)[closest]
    else
        closest = partialsortperm(norm.(distances), 1:3)
        return Vec3{Int}(sort(notes(tonnetz)[closest]))
    end
end

function trajectory_to_chord_progression(tonnetz, θ₁, θ₂, t; notemode=false)
    chords = []
    durations = []
    all_chords = project_onto.(Ref(tonnetz), θ₁, θ₂, notemode)
    for i in eachindex(all_chords)
        if i == 1 || all_chords[i] != chords[end]
            push!(chords, all_chords[i])
            if i < length(all_chords)
                push!(durations, t[i+1] - t[i])
            else
                push!(durations, t[i] - t[i-1])
            end
        else
            if i < length(all_chords)
                durations[end] += t[i+1] - t[i]
            else
                durations[end] += t[i] - t[i-1]
            end
        end
    end
    return chords, durations
end
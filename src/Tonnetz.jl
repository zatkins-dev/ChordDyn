export Tonnetz, lattice, notes, cartesian_points, angular_points

## Tonnetz struct
struct Tonnetz
    t::TriangularLattice
    notes::Matrix{Int}
end

function Tonnetz(t::TriangularLattice; startoctave::Int=3)
    notes = Matrix{Int}(undef, t.num_height, t.num_width)
    notes[:, 1] = 12startoctave .+ [0, 7, 2, -3, 4, -1, 6, 1, -4, 3, -2, 5]
    notes[:, 2] = 12startoctave .+ [4, -1, 6, 1, -4, 3, -2, 5, 0, 7, 2, -3]
    for i in 3:t.num_width
        if i % 2 == 1
            notes[:, i] = notes[:, 1] .+ (((i - 1) ÷ 2) % 12)
        else
            notes[:, i] = notes[:, 2] .+ (((i - 1) ÷ 2) % 12)
        end
        notes[:, i] -= 12 * ((notes[:, i] .+ 3) .÷ 12 .> startoctave)
    end
    return Tonnetz(t, notes)
end

function Tonnetz(; kwargs...)
    t = TriangularLattice(24, 12)
    Tonnetz(t; kwargs...)
end

lattice(t::Tonnetz) = t.t
xlims(t::Tonnetz) = xlims(lattice(t))
ylims(t::Tonnetz) = ylims(lattice(t))
xlength(t::Tonnetz) = xlength(lattice(t))
ylength(t::Tonnetz) = ylength(lattice(t))
points(t::Tonnetz) = points(lattice(t))

function cartesian_points(tonnetz::Tonnetz)
    px = [p[1] for p in points(tonnetz)[:]]
    py = [p[2] for p in points(tonnetz)[:]]
    px, py
end

function angular_points(tonnetz::Tonnetz)
    px, py = cartesian_points(tonnetz)
    bounds_x = xlims(tonnetz)
    bounds_y = ylims(tonnetz)
    2π * px / bounds_x[2], 2π * py / bounds_y[2]
end

notes(tonnetz::Tonnetz) = [n for n in tonnetz.notes[:]]
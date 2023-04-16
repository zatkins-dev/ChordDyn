export Tonnetz, lattice, notes, cartesian_points, angular_points

## Tonnetz struct
struct Tonnetz
    t::TriangularLattice
    notes::Matrix{Int}
end

function Tonnetz(t::TriangularLattice)
    notes = Matrix{Int}(undef, t.num_height, t.num_width)
    octave = [round(Int, (j - 1) / 3) + 2 for j in 1:t.num_height]
    notes[:, 1] = [(7 * (j - 1)) % 12 + 12octave[j] for j in range(1, t.num_height)]
    for i in range(2, t.num_width)
        if i % 2 == 1
            notes[:, i] = (notes[:, i-1] .- 3 .+ 12) .% 12 + 12octave
        else
            notes[:, i] = (notes[:, i-1] .+ 4) .% 12 + 12octave
        end
    end
    return Tonnetz(t, notes)
end

function Tonnetz()
    t = TriangularLattice(24, 12)
    Tonnetz(t)
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
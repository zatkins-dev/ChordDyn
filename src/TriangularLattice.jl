export TrianglularLattice, points, periodic_distance

struct TriangularLattice
    num_width::Int
    num_height::Int
    points::Matrix{Vec2f}
    xlims::Vec2f
    ylims::Vec2f
end

function TriangularLattice(num_width, num_height)
    h = 1 / num_height
    ℓ = √3 / 2 * h
    points = Matrix{Vec2}(undef, num_height, num_width)
    for i in range(1, num_width)
        points[:, i] = [Vec2f((i - 1) * ℓ, (j - 1) * h + ((i - 1) % 2) * h / 2) for j in range(1, num_height)]
    end
    return TriangularLattice(num_width, num_height, points, Vec2f(0, ℓ * num_width), Vec2f(0, 1))
end

## Getters
points(t::TriangularLattice) = t.points[:]
tri_x(t::TriangularLattice) = √3 / 2 * tri_y(t)
tri_y(t::TriangularLattice) = 1 / t.num_height
xlims(t::TriangularLattice) = t.xlims
ylims(t::TriangularLattice) = t.ylims
xlength(t::TriangularLattice) = t.xlims[2] - t.xlims[1]
ylength(t::TriangularLattice) = t.ylims[2] - t.ylims[1]

## 2D Utilities
rotation_matrix_2d(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
rotate2d(p, θ) = rotation_matrix_2d(θ) * p

rescale_to_lims(x, lims) = x * (lims[2] - lims[1]) + lims[1]

function periodic_distance(p1::Vec2f, p2::Vec2f; dims=[1, 1])
    deltas = abs.(p2 - p1)
    deltas = [min(δ, d - δ) for (δ, d) in zip(deltas, dims)]
    norm(deltas)
end



using Ordinary, ChordDyn
using DelimitedFiles
using CairoMakie
let
    t = Tonnetz()
    fig, ax = ChordDyn.plot(t)
    save("output/tonnetz.svg", fig)

    # return AssertionError("fix tonnetz plotting, especially the triangles on the edges")


    # domain = Triangle((0, 0), (xlength(tonnetz), ylength(tonnetz)))
    # recmesh = GeometryBasics.mesh(domain)
    # Makie.mesh!(recmesh; color=(:orange, 0.1))
end
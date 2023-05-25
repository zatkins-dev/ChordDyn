using GLMakie
using Makie.GeometryBasics

let
    f = Figure()
    Axis(f[1, 1])

    poly!(Polygon([Point2f(0, 0), Point2f(1, 0), Point2f(0, 1)]), color=:pink, strokecolor=:black, strokewidth=1)

    display(f)
end
using GLMakie
using EllipticFunctions
using MathTeXEngine
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/ParametricCurveTools.jl")
include("../../Tools/isothermicCylinderTools.jl")

fig = Figure(size=(1200, 800))

ax = Axis(
    fig[1, 1],
    limits =(-3, 3, -4, 2),
    title = "Isothermic Cylinder Curvature Lines - Changing Family Parameter"
)
hidedecorations!(ax)
hidespines!(ax)


animLength = 10
framerate = 30
nFrames = animLength * framerate



tau = 0.5 + 0.3 * im
omega = findOmegaRhombic(tau)
interval = LinRange(0, 2 * pi, 1000)


wIterator = LinRange(0, animLength, nFrames)

record(fig, "isothermicCylinderCurvatureLines_ParameterChanging.mp4", wIterator; framerate = framerate) do wIt
    empty!(ax)

    w = 2 * sin(wIt/3)
    curvePoints = [(rhombicLatticeCurve(t, w, omega, tau)...,) for t in interval]


    text!(
        (0.05, 0.95), 
        text= string("w = ", round(w, digits=2)), 
        fontsize = 20,
        space = :relative
    )

    text!(
        (0.05, 0.9), 
        text= "Omega = " * string(round(omega, digits=2)),
        fontsize = 20,
        space = :relative
    )

    text!(
        (0.05, 0.85), 
        text= "tau = " * string(round(real(tau), digits=2)) * " + " * string(round(imag(tau), digits=2)) * "i",
        fontsize = 20,
        space = :relative
    )
    lines!(ax, curvePoints; color=:dodgerblue, linewidth=2)
end


#=
ax = Axis(
    fig[1, 1],
    limits =(-1, 7, -4, 4),
    title = "Isothermic Cylinder Curvature Lines - Omega Parameter"
)
hidedecorations!(ax)
hidespines!(ax)


animLength = 10
framerate = 30
nFrames = animLength * framerate

tau = 0.5 + 0.3 * im
omega = findOmegaRhombic(tau)
interval = LinRange(-6 * pi, 6 * pi, 1000)
wIterator = LinRange(0.2, 0.9, nFrames)

record(fig, "isothermicCylinderCurvatureLines_OmegaChanging.mp4", wIterator; framerate = framerate) do wIt
    empty!(ax)

    w = 0.8
    curvePoints = [(rhombicLatticeCurve(t, w, wIt, tau)...,) for t in interval]

    text!(
        (0.05, 0.95), 
        text= string("w = ", round(w, digits=2)), 
        fontsize = 20,
        space = :relative
    )

    text!(
        (0.05, 0.9), 
        text= "Omega = " * string(round(wIt, digits=2)),
        fontsize = 20,
        space = :relative
    )

    text!(
        (0.05, 0.85), 
        text= "tau = " * string(round(real(tau), digits=2)) * " + " * string(round(imag(tau), digits=2)) * "i",
        fontsize = 20,
        space = :relative
    )
    lines!(ax, curvePoints; color=:dodgerblue, linewidth=2)
end
=#
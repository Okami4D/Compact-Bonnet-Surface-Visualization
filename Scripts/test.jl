using GLMakie
using EllipticFunctions
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/ParametricCurveTools.jl")
include("../Tools/isothermicCylinderTools.jl")

# Plot Setup
fig = Figure(
    size=(1200, 800), 
    scenekw = (
        lights = [DirectionalLight(RGBf(1, 1, 1), 
        Vec3f(-1, 0, 0))],
        )
)

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)
    q

tau = 0.5 + 0.3 * im
omega = findOmegaRhombic(tau)

A = 1.44531765156
B = 1.33527652772
C = 1.05005399924

w = (v) -> C + (A/pi) * sin(v) - (A/(pi^2)) *cos(v) - (B/(2 * pi)) * sin(2 * v) + (B/(4 * pi^2)) * cos(2 * v)


# theta = LinRange(0, 2* pi, 100)
# x = [cos(t) for t in theta]
# y = [sin(t) for t in theta]
# z = [w(v) for v in theta]

# lines!(x, y, z)
# fig

axisCalc = (x) -> rhombicAxisCalculation(x, omega, tau)

wFunc = (v) -> 0
rhombicAxisCalculation(0.2, omega, tau)
numericallySolveRotation(0, wFunc, axisCalc)
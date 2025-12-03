using GLMakie
import Quaternions as Q
include("../../Tools/ParametricCurveTools.jl")
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/isothermicCylinderTools.jl")

fig = Figure(size=(1200, 800))

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)

tau = 0.5 + 0.3205128205 * im # Calculated from the Paper
omega = findOmegaRhombic(tau)

A = 1.44531765156
B = 1.33527652772
C = 1.05005399924



animLength = 2
framerate = 10
nFrames = animLength * framerate

xmin =0 
xmax = 2 * pi
ymin = 0
ymax = 4 * pi

Xdensity = 200
Ydensity = 200

x = LinRange(xmin, xmax, 900)
y = LinRange(ymin, ymax, 900)

iterator = LinRange(0, ymax, nFrames)

w = (v) -> C + (A / pi) * sin(v) - (A / (pi^2)) * cos(v) - (B / (2 * pi)) * sin(2 * v) + (B / (4 * pi^2)) * cos(2 * v)

record(fig, "Out.mp4", enumerate(iterator); framerate = framerate) do (i, it)
    empty!(ax)    
    
    axis = (v) -> rhombicAxisCalculation(v, omega, tau)
    currentRotQuat = numericallySolveRotation(w, axis)(it)
    currentRot = Q.imag_part(currentRotQuat)
    f = isothermicCylinder(w, axis, omega, tau)


    arrows3d!(ax, [(0, 0, 0)], [currentRot])
    #plotParametricSurface(f, LinRange(xmin, it, Ydensity), y)
    #plotCoordinateCurves(f, x, y, [it], [])
    
    plotParametricSurface(f, x, LinRange(ymin, it, Ydensity))
    plotCoordinateCurves(f, x, y, [], [it])
    print("Frame: ", i, "/", nFrames,  "\n")
end
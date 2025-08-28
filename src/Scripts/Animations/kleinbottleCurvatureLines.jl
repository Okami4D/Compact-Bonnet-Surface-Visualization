using GLMakie

include("../../Tools/ParametricCurveTools.jl")
include("../../Tools/ParametricSurfaceTools.jl")

fig = Figure(size=(1200, 800))

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)


animLength = 4
framerate = 30
nFrames = animLength * framerate

xmin =-pi /2 
xmax = pi / 4
ymin = - pi
ymax = 2* pi / 3

Xdensity = 15
Ydensity = 20

x = LinRange(xmin, xmax, 300)
y = LinRange(ymin, ymax, 300)

iterator = LinRange(0, xmax, nFrames)


record(fig, "wenteCurvatureLines.mp4", enumerate(iterator); framerate = framerate) do (i, it)
    empty!(ax)    
    f = parametricFuncWente(12.7898)
    plotParametricSurface(f, x, LinRange(ymin, it, Ydensity))
    plotCoordinateCurves(f, x, y, [], [it])
    print("Frame: ", i, "/", nFrames,  "\n")
end

#= record(fig, "kleinbottleCurvatureLines.mp4", iterator; framerate = framerate) do it
    empty!(ax)    
    f = parametricFuncKleinBottle()
    plotParametricWireframe(f, x, y; color = (:lightgray, 0.5))
    plotCoordinateCurves(f, x, y, [it], [], xColor = :teal)
end =#
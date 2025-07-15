using GLMakie

include("../Tools/ParametricCurveTools.jl")
include("../Tools/ParametricSurfaceTools.jl")

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

f = parametricFuncKleinBottle()

xmin =0 
xmax = 2 * pi
ymin = 0
ymax = 2 * pi


x = LinRange(xmin, xmax, 1000)
y = LinRange(ymin, ymax, 1000)

density = Observable(5)  # initial density

slider = Slider(fig[2, 1], range = 0:20, startvalue = 5)

on(slider.value) do val
    density[] = val
end

on(density) do _
    empty!(ax)
    plotCoordinateCurves(f, x, y, LinRange(xmin, xmax, density[]), LinRange(ymin, ymax, density[]))
end

plotCoordinateCurves(f, x, y, LinRange(xmin, xmax, density[]), LinRange(ymin, ymax, density[]))

fig
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
    aspect = :equal, 
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


x = LinRange(xmin, xmax, 100)
y = LinRange(ymin, ymax, 100)

# initial density
e = Observable(0.0)

slider = Slider(fig[2, 1], range = 0:0.01:20, startvalue = 0)
plotParametricSurface(f, x, LinRange(ymin, e[], 10))


on(slider.value) do val
    plotCoordinateCurves(f, x, y, (val), ())
    e[] = val
end


plotParametricWireframe(f, x, y, ax; color = (:black, 0.05), transparency = true)

fig
using GLMakie
using LaTeXStrings
include("../../../Tools/ParametricSurfaceTools.jl")
include("../../../Tools/QuaternionicGeometryToolkit.jl")

# Plot Setup
fig = Figure(size=(1200, 800), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))

ax = Axis3(
    fig[1, 1], 
    aspect = :equal, 
    perspectiveness = 0.6, 
    clip=false,
    title= L"\lambda = (u, v) = 1 + 0.05 \sqrt{u^2 + v^2} \mathbb{i}"
)

hidedecorations!(ax)
hidespines!(ax)


f = parametricFuncPlane()

x = 1.0
arrow = (0.0, 1.0, 1.0)


lamda = (u, v) -> Quaternions.Quaternion(1, 0.05 * sqrt(u^2 + v^2), 0,0)


g = spinTransform(f, lamda)


# Rendering
res1 = 100
res2 = 100

u1 = LinRange(-2, 2, res1)
v1 = LinRange(-2, 2, res1)

res2 = 10
u2 = LinRange(-2, 2, res2)
v2 = LinRange(-2, 2, res2)

plotParametricSurface(f, u1, v1)
plotParametricWireframe(g, u2, v2, color = (:black, 0.1), transparency = true)

#arrows3d!( ax, (0.0, 0.0, 0.0), [arrow]; color = (:black, 0.2), normalize = true, lengthscale = 2, transparency = true)
fig
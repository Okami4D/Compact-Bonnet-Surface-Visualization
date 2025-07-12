using GLMakie
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/QuaternionicGeometryToolkit.jl")

# Plot Setup
fig = Figure(size=(1200, 800), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)


f = parametricFuncPlane()

A = 1

fu, fv = numericDifferentials(f)
df1, df2 = generateDualSurface(fu, fv)


lamda1 = (u, v) -> Quaternions.Quaternion(0, A + f(u, v)[1], A + f(u, v)[2], A + f(u, v)[3])
lamda2 = (u, v) -> Quaternions.Quaternion(0, -A + f(u, v)[1], -A + f(u, v)[2], -A + f(u, v)[3])


gu1, gv1 = spinTransform(df1, df2, lamda1)
gu2, gv2 = spinTransform(df1, df2, lamda2)

g1 = integrateForm(gu1, gv1, 0.0, 0.0)
g2 = integrateForm(gu2, gv2, 0.0, 0.0)


res1 = 100

u1 = LinRange(-2, 2, res1)
v1 = LinRange(-2, 2, res1)

res2 = 10
u2 = LinRange(-2, 2, res2)
v2 = LinRange(-2, 2, res2)


plotParametricSurface(f, u1, v1)
plotParametricWireframe(g1, u2, v2)
plotParametricWireframe(g2, u2, v2)
fig
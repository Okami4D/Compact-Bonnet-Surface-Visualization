using GLMakie
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/QuaternionicGeometryToolkit.jl")

# Plot Setup
fig = Figure(size=(1200, 800), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :equal, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)


f = parametricFuncIsothermicCatenoid()

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

u1 = LinRange(0, 2* pi, res1)
v1 = LinRange(-pi/2 + 0.01, pi/2 - 0.01, res1)

res2 = 10
u2 = LinRange(0, 2* pi, res2)
v2 = LinRange(-pi/2 + 0.01, pi/2 - 0.01, res2)


plotParametricSurface(f, u1, v1; transparency = true)
plotParametricWireframe(g1, u2, v2; color = (:red, 0.05), transparency = true)
plotParametricWireframe(g2, u2, v2; color = (:blue, 0.05), transparency = true)
fig
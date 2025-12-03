#=
This Scene visualises the isothermic Torus with one family of curvatures lines given
by the paper by Bobenko and Hoffmann. We use "bobenkoCurves" and rotate this family to
generate a surface. The rotation is calculated by the function "numericallySolveRotation".

The specific example uses a tau-admissible reparametrization calculated numerically in the paper.
=#


using GLMakie

include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/ParametricCurveTools.jl")
include("../../Tools/isothermicCylinderTools.jl")
include("../../Tools/QuaternionicGeometryToolkit.jl")


# Plot Setup
fig = Figure(
    size=(1200, 800),
    scenekw=(
        lights=[DirectionalLight(RGBf(1, 1, 1),
            Vec3f(-1, 0, 0))],
    )
)

ax = Axis3(
    fig[1, 1],
    aspect=:equal,
    perspectiveness=0.7,
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)

#-------
a = 2
#gamma = (t) -> (a * cosh(t/a), a * sinh(t/a))
gamma = (t) -> (2 + sin(5*t), t)
fps, maxS = parametricFuncIsothermicSurfaceofRevolution(gamma, (-pi, pi), nSamples = 200)
fp = parametricFuncEnneper()
A = 1

fu, fv = numericDifferentials(fp)
df1, df2 = generateDualSurface(fu, fv)


lamda1 = (u, v) -> Quaternions.Quaternion(0, A + fp(u, v)[1], A + fp(u, v)[2], A + fp(u, v)[3])
lamda2 = (u, v) -> Quaternions.Quaternion(0, -A + fp(u, v)[1], -A + fp(u, v)[2], -A + fp(u, v)[3])


gu1, gv1 = spinTransform(df1, df2, lamda1)
gu2, gv2 = spinTransform(df1, df2, lamda2)

g1 = integrateForm(gu1, gv1, 0.0, 0.0)
g2 = integrateForm(gu2, gv2, 0.0, 0.0)



umin = -2
umax = 2
vmin = -2
vmax = 2



# Basic UI
slider = Slider(fig[2, 1], range = LinRange(vmin, vmax,100), startvalue = 0)
#-------

N = 200

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)

plotParametricWireframe(g1, x, y; color = (:black, 0.05), transparency = true)
plotParametricWireframe(g2, x, y; color = (:black, 0.05), transparency = true)

curvObs = lift(slider.value) do val
    points = [f(x_val, val) for x_val in x]
end

lines!(ax, curvObs; color = (:blue, 1), linewidth = 4)

fig
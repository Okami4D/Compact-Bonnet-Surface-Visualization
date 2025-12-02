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

# Constants

tau = 0.5 + 0.3205128205 * im # Calculated from the Paper
omega = findOmegaRhombic(tau)

umin = 0 
umax = 2 * pi
vmin = 0
vmax = 4 * pi

# Magic Constants from the Paper
A = 1.44531765156
B = 1.33527652772
C = 1.05005399924
#A = 0.1
#B = 0.1
#C = 0.3



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
    perspectiveness=0.6,
    clip=false
)

ax2 = Axis3(
    fig[1, 2],
    aspect=:equal,
    perspectiveness=0.6,
    clip=false
)
ax3 = Axis3(
    fig[1, 3],
    aspect=:equal,
    perspectiveness=0.6,
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)
hidedecorations!(ax2)
hidespines!(ax2)
hidedecorations!(ax3)
hidespines!(ax3)

# Basic UI
slider = Slider(fig[2, 1], range = LinRange(vmin, vmax,100), startvalue = 0)


#-------


# Defining the reparemitriazation
w = (v) -> C + (A / pi) * sin(v) - (A / (pi^2)) * cos(v) - (B / (2 * pi)) * sin(2 * v) + (B / (4 * pi^2)) * cos(2 * v)

# Defining the axis calculation
axis = (v) -> rhombicAxisCalculation(v, omega, tau)

f = isothermicCylinder(w, axis, omega, tau)


A = 1

fu, fv = numericDifferentials(f)
df1, df2 = generateDualSurface(fu, fv)


lamda1 = (u, v) -> Quaternions.Quaternion(0, A + f(u, v)[1], A + f(u, v)[2], A + f(u, v)[3])
lamda2 = (u, v) -> Quaternions.Quaternion(0, -A + f(u, v)[1], -A + f(u, v)[2], -A + f(u, v)[3])


gu1, gv1 = spinTransform(df1, df2, lamda1)
gu2, gv2 = spinTransform(df1, df2, lamda2)

g1 = integrateForm(gu1, gv1, 0.0, 0.0)
g2 = integrateForm(gu2, gv2, 0.0, 0.0)

# Plotting the surface
N = 200

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)

plotParametricWireframe(f, x, y, ax; color = (:black, 0.05), transparency = true)

curvObs = lift(slider.value) do val
    points = [f(x_val, val) for x_val in x]
end

lines!(ax, curvObs; color = (:blue, 1), linewidth = 4)

plotParametricWireframe(g1, x, y, ax2; color = (:black, 0.05), transparency = true)
curvObs_1 = lift(slider.value) do val
    points = [g1(x_val, val) for x_val in x]
end
lines!(ax2, curvObs_1; color = (:blue, 1), linewidth = 4)
plotParametricWireframe(g2, x, y, ax3; color = (:black, 0.05), transparency = true)
curvObs_2 = lift(slider.value) do val
    points = [g2(x_val, val) for x_val in x]
end
lines!(ax3, curvObs_2; color = (:blue, 1), linewidth = 4)

fig
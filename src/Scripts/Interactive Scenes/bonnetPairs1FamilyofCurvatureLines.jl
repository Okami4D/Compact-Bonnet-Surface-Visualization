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
include("../../Tools/bonnetCalculations.jl")

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
    aspect=:data,
    perspectiveness=0.6,
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)

# Basic UI
#slider = Slider(fig[2, 1], range = LinRange(vmin, vmax,100), startvalue = 0)


#-------


# Defining the reparemitriazation
w = (v) -> C + (A / pi) * sin(v) - (A / (pi^2)) * cos(v) - (B / (2 * pi)) * sin(2 * v) + (B / (4 * pi^2)) * cos(2 * v)


f_1, f_2 = generateBonnetPair(w, omega, tau);


# Extra calculations for calculating the necessary rotation and axis to use the rotational symmetry of the isothermic cylinder
#=
rot = numericallySolveRotation(w, axis)

q = 1/rot(0) * rot(2 * pi)
rotationAngle = acos(q.s) * 2
rotationAxis =  (q.v1 / sin(rotationAngle / 2), q.v2 / sin(rotationAngle / 2), q.v3 / sin(rotationAngle / 2))
=#


# Plotting the surface
N = 100

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)

plotParametricWireframe(f_1, x, y; color = (:blue, 0.05), transparency = true)
#plotParametricWireframe(f_2, x, y; color = (:green, 0.05), transparency = true)

fig
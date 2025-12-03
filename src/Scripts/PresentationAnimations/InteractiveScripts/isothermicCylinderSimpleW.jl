#=
This Scene visualises the isothermic Torus with one family of curvatures lines given
by the paper by Bobenko and Hoffmann. We use "bobenkoCurves" and rotate this family to
generate a surface. The rotation is calculated by the function "numericallySolveRotation".

The specific example uses a tau-admissible reparametrization calculated numerically in the paper.
=#


using GLMakie
import Quaternions as Q

include("../../../Tools/ParametricSurfaceTools.jl")
include("../../../Tools/ParametricCurveTools.jl")
include("../../../Tools/isothermicCylinderTools.jl")

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

# Plot Setup
fig = Figure(
    size=(1200, 1200),
    scenekw=(
        lights=[DirectionalLight(RGBf(1, 1, 1),
            Vec3f(-1, 0, 0))],
    )
)

ax = Axis3(
    fig[2:3, 1:2],
    aspect=:data,
    perspectiveness=0.6,
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)


plotAx = Axis(
    fig[1, 1],
    limits =(0, 2 * pi, 0, 2 * pi),
)

curveAx = Axis(
    fig[1, 2],
    limits =(-2, 4, -4, 2),
)

# Basic UI

sliderGrid = GridLayout(fig[4, 1:2], rows = 3, columns = 2)
Label(sliderGrid[1, 1], "u-Coordinate =", tellwidth = false)
slider = Slider(sliderGrid[1, 2], range = LinRange(vmin, vmax,100), startvalue = 0)


#-------


# Defining the reparemitriazation
w = (v) -> 0.5 * sin(v) + 1

# Defining the axis calculation
axis = (v) -> rhombicAxisCalculation(v, omega, tau)


f = isothermicCylinder(w, axis, omega, tau)




# Plotting the surface
N = 200

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)

plotParametricWireframe(f, x, y, ax; color = (:black, 0.05), transparency = true)


curvObs = lift(slider.value) do val
    points = [f(x_val, val) for x_val in x]
end

planarCurveObs = lift(slider.value) do val
    points = [rhombicLatticeCurve(u_val, w(val), omega, tau) for u_val in x]
end

posistionCurveObs = lift(slider.value) do val
    points = [(val, 0), (val, w(val))]
end

lines!(ax, curvObs; color = (:blue, 1), linewidth = 4)

lines!(plotAx, posistionCurveObs; color = :red)
lines!(plotAx, [(t, w(t)) for t in x])

lines!(curveAx, planarCurveObs)


# Extra calculations for calculating the necessary rotation and axis to use the rotational symmetry of the isothermic cylinder

rot = numericallySolveRotation(w, axis)

q = 1/rot(0) * rot(2 * pi)
rotationAngle = acos(q.s) * 2
rotationAxis =  (q.v1 / sin(rotationAngle / 2), q.v2 / sin(rotationAngle / 2), q.v3 / sin(rotationAngle / 2))

arrows3d!(ax, [(0, 0, 0)], [rotationAxis]; color = (:black, 0.5), normalize = true, tipradius = 0.01, lengthscale = 5,shaftradius = 0.01)


fig
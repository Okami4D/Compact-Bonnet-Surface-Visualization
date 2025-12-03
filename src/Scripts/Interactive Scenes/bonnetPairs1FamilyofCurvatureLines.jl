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


# Plot Setup
fig = Figure(
    size=(1200, 800),
    scenekw=(
        lights=[DirectionalLight(RGBf(1, 1, 1),
            Vec3f(-1, 0, 0))],
    )
)

axLeft = Axis3(
    fig[1, 1],
    aspect=:equal,
    perspectiveness=0.6,
    clip=false
)


axRight = Axis3(
    fig[1, 2],
    aspect=:equal,
    perspectiveness=0.6,
    clip=false
)

hidedecorations!(axLeft)
hidedecorations!(axRight)

hidespines!(axLeft)
hidespines!(axRight)


# Basic UI
slider = Slider(fig[2, 1], range = LinRange(vmin, vmax,100), startvalue = 0)


#-------


# Defining the reparemitriazation
#w = (v) -> C + (A / pi) * sin(v) - (A / (pi^2)) * cos(v) - (B / (2 * pi)) * sin(2 * v) + (B / (4 * pi^2)) * cos(2 * v)
w = (v) -> 0.5 * sin(v) + 1
f_1, f_2 = generateBonnetPair(w, omega, tau);

# Plotting the surface
N = 80

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)

plotParametricWireframe(f_1, x, y, axLeft; color = (:black, 0.05), transparency = true)
plotParametricWireframe(f_2, x, y, axRight; color = (:black, 0.05), transparency = true)
fig


curvObs_1 = lift(slider.value) do val
    points = [f_1(x_val, val) for x_val in x]
end

curvObs_2 = lift(slider.value) do val
    points = [f_2(x_val, val) for x_val in x]
end

lines!(axLeft, curvObs_1; color = (:blue, 1), linewidth = 4)
lines!(axRight, curvObs_2; color = (:blue, 1), linewidth = 4)

fig

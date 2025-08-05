using GLMakie
using DifferentialEquations
using Quaternions

include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/ParametricCurveTools.jl")
include("../Tools/isothermicCylinderTools.jl")

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
hidedecorations!(ax)
hidespines!(ax)


tau = 0.5 + 0.3205128205 * im
omega = findOmegaRhombic(tau)

A = 1.44531765156
B = 1.33527652772
C = 1.05005399924

w = (v) -> C + (A / pi) * sin(v) - (A / (pi^2)) * cos(v) - (B / (2 * pi)) * sin(2 * v) + (B / (4 * pi^2)) * cos(2 * v)
axis = (v) -> rhombicAxisCalculation(v, omega, tau)

f = isothermicCylinder(w, axis, omega, tau)

rot = numericallySolveRotation(w, axis)

q = 1/rot(0) * rot(2 * pi)
rotationAngle = acos(q.s) * 2
rotationAxis =  (q.v1 / sin(a / 2), q.v2 / sin(a / 2), q.v3 / sin(a / 2))


function rotateFunction(f, angle, axisVector)
    axis = Quaternion(0, axisVector...)
    rotQuat = Quaternion(cos(angle/2), 0, 0, 0) + sin(angle/2) * (1/sqrt(axis.s^2 + axis.v1^2 + axis.v2^2 + axis.v3^2)) * axis
    fRotatedQuaternion = (u, v) -> rotQuat * Quaternion(0, f(u, v)...) * (1/rotQuat)
    function output(u, v)
        fRotated = fRotatedQuaternion(u, v)
        return (fRotated.v1, fRotated.v2, fRotated.v3)
    end
    return output
end


fRotated = rotateFunction(f, rotationAngle, rotationAxis)

umin = 0 
umax = 2 * pi
vmin = 0
vmax = 4 * pi

N = 200

x = LinRange(umin, umax, N)
y = LinRange(vmin, vmax, N)


plotParametricWireframe(f, x, y; color = (:black, 0.01), transparency = true)
#plotParametricWireframe(fRotated, x, y; color = (:black, 0.01), transparency = true)
slider = Slider(fig[2, 1], range = LinRange(vmin, vmax,100), startvalue = 0)

curvObs = lift(slider.value) do val
    points = [f(x_val, val) for x_val in x]
end

curvObs2 = lift(slider.value) do val
    points2 = [fRotated(x_val, val) for x_val in x]
end

lines!(ax, curvObs; color = (:blue, 1), linewidth = 4)
#lines!(ax, curvObs2; color = (:blue, 1), linewidth = 4)


arrows3d!([(0,0,0)], [rotationAxis])

fig





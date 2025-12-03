using GLMakie
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/QuaternionicGeometryToolkit.jl")
include("../../Tools/ParametricCurveTools.jl")

fig = Figure(size=(1200, 800))

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.8, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)

# Define surface to spin
f = parametricFuncEnneper()

# Parameter domain
umin, umax = -2.0, 2.0
vmin, vmax = -2.0, 2.0
udensity, vdensity = 100, 100

u = LinRange(umin, umax, udensity)
v = LinRange(vmin, vmax, vdensity)

# Animation settings
animLength = 5
framerate = 30
nFrames = animLength * framerate

# Rotation parameter: angle ∈ [0, 2π]
iterator = LinRange(0, 2π, nFrames)

record(fig, "SpinningObject.mp4", enumerate(iterator); framerate = framerate) do (i, angle)
    empty!(ax)
    angle *= 0.3
    # Rotate surface around z-axis
    f_rotated = (u_val, v_val) -> begin
        p = f(u_val, v_val)
        x, y, z = p[1], p[2], p[3]
        x_rot = x * cos(angle) - y * sin(angle)
        y_rot = x * sin(angle) + y * cos(angle)
        (x_rot, y_rot, z)
    end
    
    # Plot as wireframe with transparency
    plotParametricWireframe(f_rotated, u, v; color = (:black, 0.1), transparency = true)
    
    print("Frame: ", i, "/", nFrames, "\n")
end
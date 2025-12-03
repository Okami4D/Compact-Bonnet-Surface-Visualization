using GLMakie
import Quaternions as Q
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/QuaternionicGeometryToolkit.jl")

fig = Figure(size=(1200, 800))

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.7, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)

# Define two surfaces to morph between
f1 = parametricFuncIsothermicCatenoid()
#f2 = f1
f2 = parametricFuncIsothermicHelicoid()

# Parameter domain
umin, umax = 0, 2 * pi
vmin, vmax = -pi/2 + 0.3,  pi/2 - 0.3
udensity, vdensity = 100, 100

u = LinRange(umin, umax, udensity)
v = LinRange(vmin, vmax, vdensity)

# Animation settings
animLength = 5
framerate = 30
nFrames = animLength * framerate

# Linear interpolation parameter: t ∈ [0, 1]
iterator = LinRange(0, 1, nFrames)

record(fig, "LinearMorph.mp4", enumerate(iterator); framerate = framerate) do (i, t)
    empty!(ax)
    
    # Linear morph: f(u,v,t) = (1-t)*f1(u,v) + t*f2(u,v)
    f_morphed = (u_val, v_val) -> (1 - t) .* f1(u_val, v_val) .+ t .* f2(u_val, v_val)
    
    plotParametricSurface(f_morphed, u, v)
    
    print("Frame: ", i, "/", nFrames, "\n")
end
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
f1 = parametricFuncHelicoid(1, 1)



dir = (0, 1, 1)
rotQuat = Quaternions.Quaternion(1, dir...)

f2 = (u, v) -> Q.imag_part(Q.conj(rotQuat) * Q.Quaternion(0, f1(u, v)...) * rotQuat)

# Parameter domain
umin, umax = -1.0, 1.0
vmin, vmax = 0.0, 4.0 * pi
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
    
    arrows3d!(ax, [(0, 0, 0)], [dir])
    plotParametricSurface(f_morphed, u, v)
    
    print("Frame: ", i, "/", nFrames, "\n")
end
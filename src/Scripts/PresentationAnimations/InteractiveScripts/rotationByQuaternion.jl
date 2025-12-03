using GLMakie
import Quaternions as Q
include("../../../Tools/ParametricSurfaceTools.jl")
include("../../../Tools/QuaternionicGeometryToolkit.jl")

fig = Figure(size=(1200, 900))

ax = Axis3(
    fig[1, 1], 
    aspect = :equal, 
    perspectiveness = 0.7, 
    clip=false
)
#hidedecorations!(ax)
#hidespines!(ax)

# Define two surfaces to morph between
f1 = parametricFuncHelicoid(1, 1)

dir = (0, 1, 1)
rotQuat = Quaternions.Quaternion(1, dir...)

f2 = (u, v) -> Q.imag_part(Q.inv(rotQuat) * Q.Quaternion(0, f1(u, v)...) * rotQuat)

# Parameter domain
umin, umax = -1.0, 1.0
vmin, vmax = 0.0, 4.0 * pi
udensity, vdensity = 100, 100

u = LinRange(umin, umax, udensity)
v = LinRange(vmin, vmax, vdensity)

n = 10
xlims!(ax, -n, n)
ylims!(ax, -n, n)
zlims!(ax, 0, n * 2)

# Slider for interpolation parameter: t ∈ [0, 1]
slider = Slider(fig[2, 1], range = LinRange(0, 1, 100), startvalue = 0)

on(slider.value) do t
    empty!(ax)

    # Linear morph: f(u,v,t) = (1-t)*f1(u,v) + t*f2(u,v)
    f_morphed = (u_val, v_val) -> (1 - t) .* f1(u_val, v_val) .+ t .* f2(u_val, v_val)
    
    arrows3d!(ax, [(0, 0, 0)], [dir]; color = :red, lengthscale = 2)
    plotParametricSurface(f_morphed, u, v)
end

# Initialize with t=0
t_init = slider.value[]
f_morphed = (u_val, v_val) -> (1 - t_init) .* f1(u_val, v_val) .+ t_init .* f2(u_val, v_val)
arrows3d!(ax, [(0, 0, 0)], [dir];color = :red, lengthscale = 2)
plotParametricSurface(f_morphed, u, v)

fig
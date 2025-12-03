using GLMakie
import Quaternions as Q
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/QuaternionicGeometryToolkit.jl")

fig = Figure(size=(1200, 900))

ax = Axis3(
    fig[1, 1], 
    aspect = :equal, 
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

# Pre-compute bounds over the full morph range to fix axis limits
t_samples = LinRange(0, 1, 20)
all_points = []
for t in t_samples
    f_morphed = (u_val, v_val) -> (1 - t) .* f1(u_val, v_val) .+ t .* f2(u_val, v_val)
    for u_val in u, v_val in v
        push!(all_points, f_morphed(u_val, v_val))
    end
end

x_vals = [p[1] for p in all_points]
y_vals = [p[2] for p in all_points]
z_vals = [p[3] for p in all_points]

xlims!(ax, minimum(x_vals), maximum(x_vals))
ylims!(ax, minimum(y_vals), maximum(y_vals))
zlims!(ax, minimum(z_vals), maximum(z_vals))

# Slider for interpolation parameter: t ∈ [0, 1]
slider = Slider(fig[2, 1], range = LinRange(0, 1, 100), startvalue = 0)

on(slider.value) do t
    empty!(ax)
    
    # Linear morph: f(u,v,t) = (1-t)*f1(u,v) + t*f2(u,v)
    f_morphed = (u_val, v_val) -> (1 - t) .* f1(u_val, v_val) .+ t .* f2(u_val, v_val)
    
    arrows3d!(ax, [(0, 0, 0)], [dir])
    plotParametricSurface(f_morphed, u, v)
end

# Initialize with t=0
t_init = slider.value[]
f_morphed = (u_val, v_val) -> (1 - t_init) .* f1(u_val, v_val) .+ t_init .* f2(u_val, v_val)
arrows3d!(ax, [(0, 0, 0)], [dir])
plotParametricSurface(f_morphed, u, v)

fig
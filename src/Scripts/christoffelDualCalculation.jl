using GLMakie
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/QuaternionicGeometryToolkit.jl")

# Plot Setup
fig = Figure(size=(1800, 800), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))

# 3D axis for f
ax1 = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false,
    title = "Surface f"
)
hidedecorations!(ax1)
hidespines!(ax1)

# 2D axis for UV plane
ax2 = Axis(
    fig[1, 2],
    aspect = DataAspect(),
    title = "UV Parameter Space"
)

# 3D axis for fstar
ax3 = Axis3(
    fig[1, 3], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false,
    title = "Dual Surface f*"
)
hidedecorations!(ax3)
hidespines!(ax3)

# Surface definitions
f = parametricFuncEnneper()
fstar = generateDualSurface(f)

# Resolution parameters
res1 = 100
u1 = LinRange(0.01, 2 * pi, res1)
v1 = LinRange(0.01, pi, res1)

res2 = 30
u2 = LinRange(0.01, 4 * pi, res2)
v2 = LinRange(0.01, 1 * pi, res2)

# Plot surfaces
plotParametricSurface(f, u1, v1, ax1)
plotParametricWireframe(fstar, u2, v2, ax3; color = (:black, 0.10), transparency = true)

# UV plane setup
xlims!(ax2, u1[1], u1[end])
ylims!(ax2, v1[1], v1[end])
lines!(ax2, [u1[1], u1[end], u1[end], u1[1], u1[1]], [v1[1], v1[1], v1[end], v1[end], v1[1]], color=:gray, linewidth=2)

# Observable for the draggable point
uv_point = Observable(Point2f(0.5, 0.5))

# Corresponding 3D points
point_on_f = lift(uv_point) do pt
    Point3f(f(pt[1], pt[2])...)
end

point_on_fstar = lift(uv_point) do pt
    Point3f(fstar(pt[1], pt[2])...)
end

# Plot the draggable point in UV space
scatter!(ax2, uv_point, color=:red, markersize=15, label="Active Point")

# Plot corresponding points on both 3D surfaces
scatter!(ax1, point_on_f, color=:red, markersize=10)
scatter!(ax3, point_on_fstar, color=:red, markersize=10)

# Make the UV point draggable
deregister_interaction!(ax2, :rectanglezoom)
on(events(ax2).mouseposition) do pos
    if Makie.is_mouseinside(ax2)
        uv_point[] = pos
    end
end

fig
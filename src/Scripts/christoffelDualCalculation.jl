using GLMakie
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/QuaternionicGeometryToolkit.jl")

# Plot Setup — changed lighting a bit (two directional lights)
fig = Figure(
    size=(1800, 800),
    scenekw = (lights = [
        DirectionalLight(RGBf(0.9, 0.9, 0.9), Vec3f(-1, -1, 1)),
        DirectionalLight(RGBf(0.6, 0.6, 0.6), Vec3f(1, 0.5, -0.2)),
    ],)
)

# 3D axis for f
axf = Axis3(
    fig[1, 1],
    aspect = :data,
    perspectiveness = 0.6,
    clip = false,
    title = "Surface f"
)
hidedecorations!(axf)
hidespines!(axf)

# 3D axis for fstar
axfstar = Axis3(
    fig[1, 2],
    aspect = :data,
    perspectiveness = 0.6,
    clip = false,
    title = "Dual Surface f*"
)
hidedecorations!(axfstar)
hidespines!(axfstar)

# Surface definitions

umin, umax = -pi,  pi
vmin, vmax = -1, 1
f = parametricFuncIsothermicHelicoid()

fstar = generateDualSurface(f)

# Resolution parameters for surfaces
res1 = 50
u1 = LinRange(umin, umax, res1)
v1 = LinRange(vmin, vmax, res1)

res2 = 20
u2 = LinRange(umin, umax, res2)
v2 = LinRange(vmin, vmax, res2)

save("f2.obj", createParametricMesh(f, u1, v1))
save("fDual2.obj", createParametricMesh(fstar, u1, v1))

#=
# Plot surfaces (static)
print("Plotting surfaces...\n")
plotParametricSurface(f, u1, v1, axf)
print("Surface f plotted.\n")
plotParametricWireframe(fstar, u2, v2, axfstar; color = (:black, 0.10), transparency = true)
print("Dual surface f* plotted.\n")
# UI: sliders to control origin (fast updates only update lines on f)
u_slider = Slider(fig[2, 1], range = umin:0.01:umax, startvalue = 0, width = 500)
v_slider = Slider(fig[3, 1], range = vmin:0.01:vmax, startvalue = 0, width = 500)

# Button to compute and show the corresponding lines on fstar (expensive)
compute_button = Button(fig[4, 1])
status_label = Label(fig[4, 2], tellwidth = false)

#Line geometry parameters
# keep a single delta array starting at 0 so lines originate at the exact origin point
line_length = (umax - umin) / 4
lineSamples = 64
lineDelta = LinRange(0.0, line_length, lineSamples)

# Observables to keep track of plotted line objects so we can delete them when updating
f_lines = Observable(Vector{Any}())      # stores the two plotted objects on axf
fstar_lines = Observable(Vector{Any}())  # stores the two plotted objects on axfstar


lineU = lift(u_slider.value, v_slider.value) do u0, v0
    [(u0 + delta, v0) for delta in lineDelta]
end

lineV = lift(u_slider.value, v_slider.value) do u0, v0
    [(u0, v0 + delta) for delta in lineDelta]
end

FlineU = lift(lineU) do pts
    [Point3f(f(u, v)...) for (u, v) in pts]
end
FlineV = lift(lineV) do pts
    [Point3f(f(u, v)...) for (u, v) in pts]
end

lines!(axf, FlineU; color = :red, linewidth = 2)
lines!(axf, FlineV; color = :blue, linewidth = 2)


# Button action: compute corresponding lines on fstar (expensive) and show them in one action
on(compute_button.clicks) do _
    status_label.text[] = "Computing f* lines..."
    
    FStarlineU = [Point3f(fstar(u, v)...) for (u, v) in lineU[]]
    FStarlineV = [Point3f(fstar(u, v)...) for (u, v) in lineV[]]

    lines!(axfstar, FStarlineU; color = :red, linewidth = 2)
    lines!(axfstar, FStarlineV; color = :blue, linewidth = 2)

    status_label.text[] = "f* lines shown (computed for current origin)"
end

fig
=#
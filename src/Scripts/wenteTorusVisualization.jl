using GLMakie

include("../Tools/ParametricSurfaceTools.jl")
# --- Scene Setup ---
fig = Figure(
    size=(1200, 800),
    scenekw = (
        lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],
    )
)

ax = Axis3(
    fig[1, 1],
    aspect = :data,
    perspectiveness = 0.6,
    clip = false
)
hidedecorations!(ax)
hidespines!(ax)

# --- Sliders for u and v domains ---
u_slider = IntervalSlider(
    fig[2, 1],
    range = LinRange(-π/2, π/4, 200),
    startvalues = (-π/2, π/4),
    width = 400
)
u_slider_text = lift(u_slider.interval) do int
    "u ∈ [" * string(round(int[1], digits=2)) * ", " * string(round(int[2], digits=2)) * "]"
end
Label(fig[2, 2], u_slider_text, tellwidth=false)

v_slider = IntervalSlider(
    fig[3, 1],
    range = LinRange(-π, 2π/3, 200),
    startvalues = (-π, 2π/3),
    width = 400
)
v_slider_text = lift(v_slider.interval) do int
    "v ∈ [" * string(round(int[1], digits=2)) * ", " * string(round(int[2], digits=2)) * "]"
end
Label(fig[3, 2], v_slider_text, tellwidth=false)

# --- Observable for the surface points ---
wente = parametricFuncWente(12.7898)



surface_obs = lift(u_slider.interval, v_slider.interval) do u_int, v_int
    u = LinRange(u_int..., 100)
    v = LinRange(v_int..., 100)
    [wente(u_val, v_val) for u_val in u, v_val in v]
end



surface!(ax, surface_obs[]...; color = (:lightgray, 0.5), shading = true)

fig
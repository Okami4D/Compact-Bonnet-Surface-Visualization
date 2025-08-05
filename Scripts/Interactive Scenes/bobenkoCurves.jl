using GLMakie
using EllipticFunctions
include("../../Tools/ParametricSurfaceTools.jl")
include("../../Tools/ParametricCurveTools.jl")
include("../../Tools/isothermicCylinderTools.jl")

# Constants
tau = 0.5 + 0.3 * im




# Plot Setup
fig = Figure(size=(1200, 800))

ax = Axis(
    fig[1, 1],
    limits =(-2, 4, -4, 2),
    title = "Isothermic Cylinder Curvature Lines - Changing Family Parameter"
)
hidedecorations!(ax)
hidespines!(ax)

sliderGrid = GridLayout(fig[2, 1], columns = 2, rows = 4, showaxis=false)

familySlide = Slider(
    sliderGrid[1, 2], 
    range = 0:0.01:2, 
    startvalue = 0.5, 
    update_while_dragging=true, 
)
familySlideText = lift(familySlide.value) do num
    string("w = ", round.(num, digits = 2))
end
Label(sliderGrid[1, 1], familySlideText, tellwidth = false)

intervalSlide = IntervalSlider(
    sliderGrid[2, 2], 
    range = LinRange(-5 *pi, 5 * pi, 1000), 
    startvalues = (0, 2 * pi), 
)
intervalSlideText = lift(intervalSlide.interval) do int
    string("Interval: ", round.(int, digits = 2))
end
Label(sliderGrid[2, 1], intervalSlideText, tellwidth = false)


omegaSlide = Slider(
    sliderGrid[3, 2], 
    range = 0:0.001:2, 
    startvalue = findOmegaRhombic(tau), 
    update_while_dragging=true, 
)
omegaSlideText = lift(omegaSlide.value) do int
    string("Omega = ", round.(int, digits = 2))
end
Label(sliderGrid[3, 1], omegaSlideText, tellwidth = true)

togl = Observable{Bool}(true)
rhombicToggle = ToggleButton(
    sliderGrid[4, 2], 
    label = "Rhombic Lattice", 
    startvalue = true, 
    width = 150,
)

on(rhombicToggle.value) do val
    togl[] = val
end

curve_obs = lift(familySlide.value, intervalSlide.interval, omegaSlide.value) do w_val, int, omegaVal
    t_vals = LinRange(int..., 1000)
    points = [(
        rhombicLatticeCurve(t, w_val, omegaVal, tau)...,
    ) for t in t_vals]
    points
end

lines!(ax, curve_obs; color=:black, linewidth=3)

#=
curve_obs = lift(familySlide.value, intervalSlide.interval, omegaSlide.value) do w_val, int, omegaVal
    t_vals = LinRange(int..., 1000)
    points = [(
        rectangularLatticeCurve(t, w_val, omegaVal, tau)...,
        0
    ) for t in t_vals]
    points
end

lines!(ax, curve_obs; color=:orange, linewidth=3)
=#

fig
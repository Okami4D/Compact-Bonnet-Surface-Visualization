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
    title = "Isothermic Cylinder Curvature Lines"
)
#hidedecorations!(ax)
hidespines!(ax)

# A grid layout for the sliders and labels
sliderGrid = GridLayout(fig[2, 1], rows = 4, columns = 2)

## UI Setup

# The slider changing the w parameter, which should be thought of as chosing one curve from the family of isothermic cylinder curves
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


# The slider changing the interval of the parameter t, which is used to plot the curve
intervalSlide = IntervalSlider(
    sliderGrid[2, 2], 
    range = LinRange(-5 *pi, 5 * pi, 1000), 
    startvalues = (0, 2 * pi), 
)
intervalSlideText = lift(intervalSlide.interval) do int
    string("Interval: ", round.(int, digits = 2))
end
Label(sliderGrid[2, 1], intervalSlideText, tellwidth = false)

# The slider changing the omega parameter, which is used to define the lattice
omegaSlide = Slider(
    sliderGrid[3, 2], 
    range = 0:0.001:4, 
    startvalue = findOmegaRhombic(tau), 
    update_while_dragging=true, 
)
omegaSlideText = lift(omegaSlide.value) do int
    string("Omega = ", round.(int, digits = 2))
end
Label(sliderGrid[3, 1], omegaSlideText, tellwidth = false)


# The toggle to switch between rhombic and rectangular lattice
rhombicToggle = Toggle(
    sliderGrid[4, 2], 
    active = toggl[]
)
togglText = lift(rhombicToggle.active) do int
    if int == true
        "Rhombic Lattice"
    else
        "Rectangular Lattice"
    end
end
Label(sliderGrid[4, 1], togglText, tellwidth = false)




# Updating the intitial omega value based on the toggle state so that thee curves close
on(rhombicToggle.active) do val
    if val == true
        omegaSlide.value = findOmegaRhombic(tau)
    else
        omegaSlide.value = findOmegaRectangular(tau)
    end
end


# Adding the slider grid to the figure
colsize!(sliderGrid, 2, Relative(0.8))


# A collection of points for the isothermic cylinder curve
curve_obs = lift(familySlide.value, intervalSlide.interval, omegaSlide.value, rhombicToggle.active) do w_val, int, omegaVal, rhombic
    t_vals = LinRange(int..., 1000)
    if rhombic
        points = [(
            rhombicLatticeCurve(t, w_val, omegaVal, tau)...,
        ) for t in t_vals]
    else
        points = [(
            rectangularLatticeCurve(t, w_val, omegaVal, tau)...,
        ) for t in t_vals]
    end
    points
end

lines!(ax, curve_obs; color=:black, linewidth=3)

fig
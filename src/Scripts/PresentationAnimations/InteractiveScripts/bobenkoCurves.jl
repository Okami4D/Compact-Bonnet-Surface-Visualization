using GLMakie
using EllipticFunctions
using LaTeXStrings

include("../../../Tools/ParametricSurfaceTools.jl")
include("../../../Tools/ParametricCurveTools.jl")
include("../../../Tools/isothermicCylinderTools.jl")

# Constants
# tau = 0.5 + 0.3 * im
tau_real = 0.5
tau_imag_init = 0.3

# Plot Setup
fig = Figure(size=(1200, 800))

ax = Axis(
    fig[1, 1],
    limits =(-2, 4, -4, 2),
    title = L"\gamma(u, w) = -i \frac{2 \vartheta_{2}(\omega)^{2}}{\vartheta_{1}'(0)\vartheta_{1}(2\omega)} \frac{\vartheta_{1}\left( \frac{1}{2} (u + iw - 3\omega)\right)}{\vartheta_{1}\left( \frac{1}{2} (u + iw + \omega)\right)} e^{(u + iw)} \frac{\vartheta_{2}' (\omega)}{\vartheta_{2} (\omega)}"
)
#hidedecorations!(ax)
hidespines!(ax)

# A grid layout for the sliders and labels
# increase rows to accomodate the extra tau-imag slider
sliderGrid = GridLayout(fig[2, 1], rows = 5, columns = 2)

## UI Setup

# The slider changing the w parameter, which should be thought of as chosing one curve from the family of isothermic cylinder curves
familySlide = Slider(
    sliderGrid[1, 2], 
    range = 0:0.01:2*pi, 
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
# startvalue uses initial tau_imag_init
omegaSlide = Slider(
    sliderGrid[3, 2], 
    range = 0:0.001:4, 
    startvalue = findOmegaRhombic(tau_real + tau_imag_init*im), 
    update_while_dragging=true, 
)
omegaSlideText = lift(omegaSlide.value) do int
    string("Omega = ", round.(int, digits = 2))
end
Label(sliderGrid[3, 1], omegaSlideText, tellwidth = false)


# New: slider to control the imaginary part of tau
imagTauSlide = Slider(
    sliderGrid[4, 2],
    range = 0:0.001:0.3547,
    startvalue = tau_imag_init,
    update_while_dragging=true,
)
imagTauText = lift(imagTauSlide.value) do v
    string("Im(tau) = ", round(v, digits=4))
end
Label(sliderGrid[4, 1], imagTauText, tellwidth = false)


# The toggle to switch between rhombic and rectangular lattice
rhombicToggle = Toggle(
    sliderGrid[5, 2], 
    active = true
)
togglText = lift(rhombicToggle.active) do int
    if int == true
        "Rhombic Lattice"
    else
        "Rectangular Lattice"
    end
end
Label(sliderGrid[5, 1], togglText, tellwidth = false)




# Updating the intitial omega value based on the toggle state so that thee curves close
on(rhombicToggle.active) do val
    current_tau = tau_real + imagTauSlide.value[] * im
    if val == true
        omegaSlide.value = findOmegaRhombic(current_tau)
        ax.title = L"\gamma(u, w) = -i \frac{2 \vartheta_{2}(\omega)^{2}}{\vartheta_{1}'(0)\vartheta_{1}(2\omega)} \frac{\vartheta_{1}\left( \frac{1}{2} (u + iw - 3\omega)\right)}{\vartheta_{1}\left( \frac{1}{2} (u + iw + \omega)\right)} e^{(u + iw)} \frac{\vartheta_{2}' (\omega)}{\vartheta_{2} (\omega)}"
    else
        omegaSlide.value = findOmegaRectangular(current_tau)
        ax.title = L"\gamma(u, w) = -i \frac{2 \vartheta_{4}(\omega)^{2}}{\vartheta_{1}'(0)\vartheta_{1}(2\omega)} \frac{\vartheta_{1}\left( \frac{1}{2} (u + iw - 3\omega)\right)}{\vartheta_{1}\left( \frac{1}{2} (u + iw + \omega)\right)} e^{(u + iw)} \frac{\vartheta_{4}' (\omega)}{\vartheta_{4} (\omega)}"
    end
end

# Also update omega when the imaginary part of tau changes
on(imagTauSlide.value) do imval
    current_tau = tau_real + imval * im
    if rhombicToggle.active[]
        omegaSlide.value = findOmegaRhombic(current_tau)
    else
        omegaSlide.value = findOmegaRectangular(current_tau)
    end
end


# Adding the slider grid to the figure
colsize!(sliderGrid, 2, Relative(0.8))


# A collection of points for the isothermic cylinder curve
curve_obs = lift(familySlide.value, intervalSlide.interval, omegaSlide.value, rhombicToggle.active, imagTauSlide.value) do w_val, int, omegaVal, rhombic, imagTau
    t_vals = LinRange(int..., 1000)
    tau = tau_real + imagTau * im
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
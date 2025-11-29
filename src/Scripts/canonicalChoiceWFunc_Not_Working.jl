using GLMakie
using Roots
using EllipticFunctions
using ForwardDiff
using FiniteDifferences
using ArbNumerics

include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/ParametricCurveTools.jl")
include("../Tools/isothermicCylinderTools.jl")
include("../Tools/jacboiThetaExtension.jl")

function complex_step_derivative(f, x, h = 1e-20)
    return imag(f(x + h * im)) / h
end

# Constants
tau = 0.5 + 0.3205128205 * im
omega = findOmegaRhombic(tau)

s1 = -3.13060628
s2 = 0.5655771591
delta = 1.61245155


umin = 0 
umax = 2 * pi
vmin = 0
vmax = 4 * pi

Q3 = (s) -> sphericalrhombicLinesQ3(tau)(s)
Q3prime = (s) -> sphericalrhombicLinesQ3Prime(tau)(s)
Q3prime2 = (s) -> sphericalrhombicLinesQ3PrimePrime(tau)(s)

gQ3 = ellipticCurveinvariantsQ3(tau)
hpQ3 = halfPeriods(gQ3...)
tauQ3 = hpQ3[2] / hpQ3[1]


QCurve = (s) -> sphericalrhombicLinesQ(s1, s2, delta, tau)(s)
QCurvePrime = (s) -> sphericalrhombicLinesQPrime(s1, s2, delta, tau)(s)
QCurvePrime2 = (s) -> sphericalrhombicLinesQPrimePrime(s1, s2, delta, tau)(s)

gQ = ellipticCurveInvariantsdeltaQ(s1, s2, delta, tau)
hpQ = halfPeriods(gQ...)

s1min = Roots.find_zero((s) -> real(QCurve(s)), s1 - delta)
s1plu = Roots.find_zero((s) -> real(QCurve(s)), s2 + delta)

if s1min > s1plu
    s1min, s1plu = s1plu, s1min
end

s0 = (th1(omega, tauQ3)^2)/(th2(omega, tauQ3)^2)



a0 = 1/24 * Q3prime2(s0)
a1 = 1/4 * Q3prime(s0)
a2 = -1/96 * Q3prime(s0) * delta^(-2) * QCurvePrime2(s1min)
a3 = (s1min - s0)
a4 = 1/4 * delta^(-2) * QCurvePrime(s1min) - (s1min - s0) * 1/24 * delta^(-2) * QCurvePrime2(s1min)

function w(s)
    wpval = weierstrass_p(ArbComplex(s), ArbComplex(hpQ3[1]), ArbComplex(hpQ3[2]))          # ℘(s)
    numerator   = a1 * wpval + a2
    denominator = a3 * wpval + a4
    z = ArbComplex(a0 + numerator / denominator)

    return weierstrass_inv_p(z, ArbComplex(hpQ[1]), ArbComplex(hpQ[2]))
end

dwFunc = (s) -> central_fdm(5, 1)(w, s)





fig = Figure(size=(1200, 800))
ax = Axis(fig[1, 1])

x = LinRange(vmin, vmax, 1000)



lines!(
    ax,
    x,
    [real(w(v)) for v in x],
    color = :blue
)

lines!(
    ax,
    x,
    [imag(w(v)) for v in x],
    color = :red
)

fig



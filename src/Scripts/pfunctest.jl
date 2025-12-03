using GLMakie
using ArbNumerics
using EllipticFunctions
using Elliptic
using PyCall

include("../Tools/jacboiThetaExtension.jl")
include("../Tools/ParametricSurfaceTools.jl")
include("../Tools/ParametricCurveTools.jl")
include("../Tools/isothermicCylinderTools.jl")


ArbNumerics.setprecision(64)

function complex_step_derivative(f, x, h = 1e-20)
    return imag(f(x + h * im)) / h
end

fig = Figure(size=(1200, 2400))
axTop = Axis3(
    fig[1, 1],
    aspect = :equal, 
    perspectiveness = 0.6
)

axBot = Axis3(
    fig[2, 1],
    aspect = :equal, 
    perspectiveness = 0.6
)


tau = 0.5 + 0.3205128205 * im
omega = findOmegaRhombic(tau)

s1 = -3.13060628
s2 = 0.5655771591
delta = 1.61245155



n = 3
res = 200

realVals = LinRange(-n, n, res)
imagVals = LinRange(-n, n, res)

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
tauQ = hpQ[2] / hpQ[1]

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


surface!(
    axTop,
    realVals,
    imagVals,
    (x, y) -> abs(Q3(x + im * y)),
    colormap = :amp,
    transparency = true
)


#wp(2z, omega = hpQ)

surface!(
    axBot,
    realVals,
    imagVals,
    (x, y) -> abs(QCurve(x + im * y)),
    colormap = :amp,
    transparency = true
)


fig


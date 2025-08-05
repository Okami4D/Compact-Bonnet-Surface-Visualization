import Quaternions: Quaternion

using Roots
using Integrals
using FiniteDifferences
using DifferentialEquations
using Zygote
using SciMLBase

include("../Tools/jacboiThetaExtension.jl")

function findOmegaRectangular(tau, n =10)
    fReal = (omega) -> real(th4prime(omega, tau))
    fImag = (omega) -> imag(th4prime(omega, tau))
    solReal = Roots.find_zero(fReal, 1)
    return solReal
end

function findOmegaRhombic(tau, n =10)
    fReal = (omega) -> real(th2prime(omega, tau))
    fImag = (omega) -> imag(th2prime(omega, tau))
    solReal = Roots.find_zero(fReal, 1)
    return solReal
end


function rectangularLatticeCurve(u, w, omega, tau = 0.5 + 25/78 * im) 
    coeff = (2 * (th4(omega, tau))^2) * (-1 * im)/((th1prime(0, tau)) * th1(2 * omega, tau))

    complexOut = coeff * th1(0.5 * (u + im * w - 3 * omega), tau)/th1(0.5 * (u + im * w + omega), tau)
    out = complexOut * exp((u + im * w) * (th4prime(omega, tau)/th4(omega, tau)))
    return (real(out), imag(out))
end

function rhombicLatticeCurve(u, w, omega, tau = 0.5 + 25/78 * im) 
    coeff = (2 * (th2(omega, tau))^2) * (-1 * im)/((th1prime(0, tau)) * th1(2 * omega, tau))

    complexOut = coeff * th1(0.5 * (u + im * w - 3 * omega), tau)/th1(0.5 * (u + im * w + omega), tau)
    out = complexOut * exp((u + im * w) * (th2prime(omega, tau)/th2(omega, tau)))
    return (real(out), imag(out))
end

function rectangularAxisCalculation(w, omega, tau)
    val1 = im * th1prime(0, tau)/(2 * th4(omega, tau)) * th4(omega - im * w, tau)/th1(im * w, tau)
    value = val1 * exp((im * w) * (th4prime(omega, tau)/th4(omega, tau))) 
    return Quaternion(real(value), imag(value), 0, 0)
end

function rhombicAxisCalculation(w, omega, tau)
    val1 = im * th1prime(0, tau)/(2 * th2(omega, tau)) * th2(omega - im * w, tau)/th1(im * w, tau)
    value = val1 * exp((im * w) * (th2prime(omega, tau)/th2(omega, tau)))
    return Quaternion(real(value), imag(value), 0, 0)
end

function numericallySolveRotation(wFunc, axisCalc)
    dwFunc = (v) -> something(gradient(wFunc, v)[1], 0)

    function Q_rhs(v)
        s = sqrt(1 - (dwFunc(v)^2))
        W = axisCalc(wFunc(v))

        return s * W * Quaternion(0.0, 0.0, 0.0, 1.0)  # Pure imaginary in k direction
    end

    function ODE!(du, u, p, t)
        Qf = Q_rhs(t)
        PHI = Quaternion(u[1], u[2], u[3], u[4])  # Convert to Q.Quaternion
        dPHI = Qf * PHI  # Q.Quaternion multiplication
        du[1] = dPHI.s
        du[2] = dPHI.v1
        du[3] = dPHI.v2
        du[4] = dPHI.v3
    end

    vspan = (0.0, 4 * pi)           # Solve from v = 0 to 2π

    Phi0 = [1.0, 0.0, 0.0, 0.0]  # Initial condition for the quaternion
    prob = ODEProblem(ODE!, Phi0, vspan)
    sol = solve(prob)

    outputFunction = (v) -> Quaternion(sol(v)[1], sol(v)[2], sol(v)[3], sol(v)[4])
    return outputFunction
end


function isothermicCylinder(w, axisCalc, omega = findOmegaRhombic(0.5 + 25/78 * im), tau = 0.5 + 25/78 * im)
    rotationFunc = numericallySolveRotation(w, axisCalc)

    function cylinderCalc(u, v)
        rot = rotationFunc(v)
        invRot = 1/rotationFunc(v)

        embeddCurve = Quaternion(rhombicLatticeCurve(u, w(v), omega, tau)..., 0, 0) * Quaternion(0, 0, 1, 0)

        calculation = invRot * embeddCurve * rot
        return (calculation.v1, calculation.v2, calculation.v3)
    end

    return cylinderCalc
end


function sphericalrhombicLinesQ3Coefficients(tau)
    omega = findOmegaRhombic(tau)
    coeff1 = (th1prime(0, tau)/th2(omega, tau))^2
    coeff2 = (th1(omega, tau)/ th2(0, tau))^2
    coeff3 = (th3(omega, tau)/ th4(0, tau))^2
    coeff4 = (th4(omega, tau)/ th3(0, tau))^2

    return coeff1, coeff2, coeff3, coeff4
end

function sphericalrhombicLinesQ3(tau)
    coeff1, coeff2, coeff3, coeff4 = sphericalrhombicLinesQ3Coefficients(tau)
    return (s) -> coeff1 * (s - coeff2) * (s - coeff3) * (s - coeff4)
end

function sphericalrhombicLinesQ(s1, s2, delta, tau)
    Q3 = sphericalrhombicLinesQ3(tau)
    return (s) -> -(s - s1)^2 * (s - s^2)^2 + delta^2 * Q3(s)
end


## Calcualted by Hand. Doublecheck maybe
function ellipticCurveinvariantsQ3(tau)
    A, B, C, D = sphericalrhombicLinesQ3Coefficients(tau)
    c1 = A/4
    c2 = -(C + B - A * D)/6
    c3 = (D * C + D * B + B * C) /4
    c4 = -(B * C * D)

    g2 = - 4 * c1 * c3 + 3 * c2^2
    g3 = 2 * c1 * c2 * c3 - c2^3 - c1^2 * c4 
    return g2, g3
end

## not sure how to read the paper since s1 and s2 are specified abut s0 is used which seems independent of s1 and s2
function sphericalCurvatureLinesW(s1, s2, delta, tau)
    Q3 = sphericalrhombicLinesQ3(tau)
    Q = sphericalrhombicLinesQ(s1, s2, delta, tau)

    Q3prime = (s) -> something(gradient(Q3, s)[1], 0)
    Qprime = (s) -> something(gradient(Q, s)[1], 0)

    Q3primeprime = (s) -> something(gradient(Q3prime, s)[1], 0)
    Qprimeprime = (s) -> something(gradient(Qprime, s)[1], 0)

    a0 = 1/24 * Q3primeprime(s)
end
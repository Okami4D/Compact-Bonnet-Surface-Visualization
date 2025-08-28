
using Roots
using Integrals
using FiniteDifferences
using DifferentialEquations
using Zygote
using SciMLBase
using ForwardDiff

import Quaternions as Q

include("../Tools/jacboiThetaExtension.jl")

function findOmegaRectangular(tau)
    fReal = (omega) -> real(th4prime(omega, tau))
    fImag = (omega) -> imag(th4prime(omega, tau))
    solReal = Roots.find_zero(fReal, 0.7)
    solImag = Roots.find_zero(fImag, 0.5)
    return solReal
end

function findOmegaRhombic(tau)
    fReal = (omega) -> real(th2prime(omega, tau))
    fImag = (omega) -> imag(th2prime(omega, tau))
    solReal = Roots.find_zero(fReal, 0.7)
    solImag = Roots.find_zero(fImag, 0.5)
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
    return Q.Quaternion(real(value), imag(value), 0, 0)
end

function rhombicAxisCalculation(w, omega, tau)
    val1 = im * th1prime(0, tau)/(2 * th2(omega, tau)) * th2(omega - im * w, tau)/th1(im * w, tau)
    value = val1 * exp((im * w) * (th2prime(omega, tau)/th2(omega, tau)))
    return Q.Quaternion(real(value), imag(value), 0, 0)
end

function numericallySolveRotation(wFunc, axisCalc)
    #dwFunc = (v) -> ForwardDiff.derivative(wFunc, v)
    dwFunc = (s) -> round(central_fdm(5, 1)(wFunc, s), digits = 20)^2

    function Q_rhs(v)
        s =  sqrt(1 - dwFunc(v))
        W = axisCalc(wFunc(v))

        return s * W * Q.Quaternion(0.0, 0.0, 0.0, 1.0)  # Pure imaginary in k direction
    end

    function ODE!(du, u, p, t)
        Qf = Q_rhs(t)
        PHI = Q.Quaternion(u[1], u[2], u[3], u[4])  # Convert to Quaternion
        dPHI = Qf * PHI  # Quaternion multiplication
        du[1] = dPHI.s
        du[2] = dPHI.v1
        du[3] = dPHI.v2
        du[4] = dPHI.v3
    end

    vspan = (0.0, 4 * pi)           # Solve from v = 0 to 2π

    Phi0 = [1.0, 0.0, 0.0, 0.0]  # Initial condition for the quaternion
    prob = ODEProblem(ODE!, Phi0, vspan)
    sol = solve(prob)

    outputFunction = (v) -> Q.Quaternion(sol(v)[1], sol(v)[2], sol(v)[3], sol(v)[4])
    return outputFunction
end


function isothermicCylinder(w, axisCalc, omega = findOmegaRhombic(0.5 + 25/78 * im), tau = 0.5 + 25/78 * im)
    rotationFunc = numericallySolveRotation(w, axisCalc)

    function cylinderCalc(u, v)
        rot = rotationFunc(v)
        invRot = 1/rotationFunc(v)

        embeddCurve = Q.Quaternion(rhombicLatticeCurve(u, w(v), omega, tau)..., 0, 0) * Q.Quaternion(0, 0, 1, 0)

        calculation = invRot * embeddCurve * rot
        return (calculation.v1, calculation.v2, calculation.v3)
    end

    return cylinderCalc
end

#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

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
    return (s) -> -1 * (s - s1)^2 * (s - s2)^2 + delta^2 * Q3(s)
end


## Calcualted by Hand. Doublecheck maybe
function ellipticCurveinvariantsQ3(tau)
    Q3_A, Q3_B, Q3_C, Q3_D = sphericalrhombicLinesQ3Coefficients(tau)
    
    A = -(Q3_A * Q3_B * Q3_C * Q3_D)
    B = Q3_A * Q3_C * Q3_D + Q3_A * Q3_B * Q3_D + Q3_A * Q3_B * Q3_C
    C = -(Q3_A * Q3_C + Q3_A * Q3_B + Q3_A  * Q3_D)
    D = Q3_A
    E = 0

    return calculateInvariants(A, B, C, D, E)
end

## Calculate the invariants of the elliptic curve Q(s) = A + B * s + C * s^2 + D * s^3 + E * s^4
function calculateInvariants(A, B, C, D, E)
    c0 = E
    c1 = D/4
    c2 = C/6
    c3 = B/4
    c4 = A

    g2 = c0 * c4 - 4 * c1 * c3 + 3 * c2^2
    g3 = c0 * c2 * c4 + 2 * c1 * c2 * c3 - c2^3 - c0 * c3^2 - c1^2 * c4 
    return g2, g3
end

function ellipticCurveInvariantsdeltaQ(s1, s2, delta, tau)
    Q3_A, Q3_B, Q3_C, Q3_D = sphericalrhombicLinesQ3Coefficients(tau)
    
    A = -(Q3_A * Q3_B * Q3_C * Q3_D)
    B = Q3_A * Q3_C * Q3_D + Q3_A * Q3_B * Q3_D + Q3_A * Q3_B * Q3_C
    C = -(Q3_A * Q3_C + Q3_A * Q3_B + Q3_A  * Q3_D)
    D = Q3_A
    E = 0

    invDelta = 1 / delta^2

    Ap = -invDelta * (s1^2 * s2^2)
    Bp = 2 * invDelta * (s1^2 * s2 + s1 * s2^2)
    Cp = - invDelta * (s1^2 + 4 * s1 * s2 + s2^2)
    Dp = 2 * invDelta * (s1 + s2)
    Ep = -invDelta

    return calculateInvariants(Ap + A, Bp + B, Cp + C, Dp + D, Ep + E)
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
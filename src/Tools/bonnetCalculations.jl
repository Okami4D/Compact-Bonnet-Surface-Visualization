
using Roots
using Integrals
using FiniteDifferences
using DifferentialEquations
using Zygote
using SciMLBase
using ForwardDiff

import Quaternions as Q

include("../Tools/jacboiThetaExtension.jl")
include("../Tools/isothermicCylinderTools.jl")

function generateBonnetPair(w, omega, tau, epsilon = 1)
    Rraw = (2 * th2(omega, tau))^2/(th1prime(0, tau) * th1(2 * omega, tau))
    R = Q.quat(real(Rraw), imag(Rraw), 0.0, 0.0) #Just use this

    axis = (v) -> rhombicAxisCalculation(v, omega, tau)

    function BHat(u, w)
        coefficient = Rraw^2 * (th1(2 * omega, tau))/(th1prime(0, tau))
        firstTerm = (th2primeprime(omega, tau))/(th2(omega, tau)) * w/2
        tempInput = (u + im * w - omega)/2
        secondTerm = imag(th2prime(tempInput, tau)/th2(tempInput, tau))

        output = coefficient * (firstTerm - secondTerm)
        return Q.Quaternion(real(output), imag(output), 0.0, 0.0)
    end


    Phi = numericallySolveRotation(w, axis)
    invPhi = (v) -> 1/Phi(v) #look into quaternion inverse

    BTilde = numericallySolveBTilde(Phi, R, w, tau)

    f = (u, v) -> Q.quat(0.0, isothermicCylinder(w, axis, omega, tau)(u, v)...)

    firstTermOut = (u, v) -> R^2 * f(pi - 2 * omega + u, v) - epsilon^2 * f(pi - u, v)
    secondTermOut = (u, v) -> 2 * epsilon * ((invPhi(v) * BHat(u, w(v)) * Q.quat(0.0 ,1.0, 0.0, 0.0) * Phi(v)) + BTilde(v))

    out_pos = (u, v) -> Q.imag_part(firstTermOut(u, v) + secondTermOut(u, v))
    out_neg = (u, v) -> Q.imag_part(firstTermOut(u, v) - secondTermOut(u, v))
    return (out_pos, out_neg)
end

function numericallySolveBTilde(Phi, R, wFunc, tau)
    #dwFunc = (v) -> ForwardDiff.derivative(wFunc, v)
    dwFunc = (s) -> round(central_fdm(5, 1)(wFunc, s), digits = 20)^2

    function Q_rhs(v)
        s =  sqrt(1 - dwFunc(v))
        b = btilde(wFunc(v), omega, tau)
        invPhi = (v) -> 1/Phi(v)
        return invPhi(v) * ( R * s * b * Q.Quaternion(0.0, 0.0, 0.0, 1.0)) * Phi(v)  # Pure imaginary in k direction
    end

    function ODE!(du, u, p, t)
        Qf = Q_rhs(t)
        dBHAT = Qf  # Quaternion multiplication
        du[1] = dBHAT.s
        du[2] = dBHAT.v1
        du[3] = dBHAT.v2
        du[4] = dBHAT.v3
    end

    vspan = (0.0, 4 * pi)           # Solve from v = 0 to 2π

    Phi0 = [0.0, 0.0, 0.0, 0.0]  # Initial condition for the quaternion
    prob = ODEProblem(ODE!, Phi0, vspan)
    sol = solve(prob)

    outputFunction = (v) -> Q.Quaternion(sol(v)[1], sol(v)[2], sol(v)[3], sol(v)[4])
    return outputFunction
end

function btilde(w, omega, tau)
    coefficient = (2 * th2(omega, tau))/(th1prime(0, tau)) * (th2(im * w - omega, tau))/(th1(im * w, tau))
    firstTerm = coefficient * (((th2primeprime(omega, tau))/(th2(omega, tau)) * w/2) - imag(th2prime(im * w/2, tau)/th2(im * w/2, tau)))
    secondTerm = im * (th1(im * w /2 - omega, tau))^2/(th2(im * w /2, tau))^2

    output = firstTerm - secondTerm
    R = real(output)
    I = imag(output)
    return Q.Quaternion(R, I, 0.0, 0.0)
end
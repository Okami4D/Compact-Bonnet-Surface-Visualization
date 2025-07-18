using Roots
using Integrals
using Quaternions
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



# VERY WIP
function numericallySolveRotation(v, wFunc, axisCalc)
    dwFunc = something(gradient(wFunc, v)[1], 0)
    # ========== 2. Define right-hand side quaternion Q(v) ==========
    # Q(v) = sqrt(1 - w'(v)^2) * W1(w(v)) * k
    function Q_rhs(v)
        s = sqrt(1 - (dwFunc^2))
        wval = wFunc(v)
        W = axisCalc(wval)
        return s * W * Quaternion(0.0, 0.0, 0.0, 1)  # Pure imaginary in k direction
    end

    function ODE(du, u, p, t)
        Q = Q_rhs(t)
        PHI = Quaternion(u[1], u[2], u[3], u[4])  # Convert to Quaternion
        dPHI = Q * PHI  # Quaternion multiplication
        du[1] = dPHI.s
        du[2] = dPHI.v1
        du[3] = dPHI.v2
        du[4] = dPHI.v3
    end

    Phi0 = [0.0, 0.0, 0.0, 1.0]  # Represented as vector: [s, i, j, k]
    vspan = (0.0, 2pi)           # Solve from v = 0 to 2π

    prob = ODEProblem(ODE, Phi0, vspan)
    sol = solve(prob, Tsit5())
    return Quaternion(sol(v)[1], sol(v)[2], sol(v)[3], sol(v)[4])  # Convert back to Quaternion
end


function isothermicCylinder(w, axisCalc, omega = findOmegaRhombic(0.5 + 25/78 * im), tau = 0.5 + 25/78 * im)

    wFunc = (v) -> w(v)
    sol = (v) -> numericallySolveRotation(v, wFunc, axisCalc)

    curve = (u, v) -> Quaternion(0, rhombicLatticeCurve(u, v, omega, tau)..., 0)

    f = (u, v) -> 1/sol(v) * curve(u, w(v)) * Quaternion(0.0, 0.0, 1.0, 0.0) * sol(v)

    # 4. Calculate the coordinates of the isothermic cylinder
    return f
end
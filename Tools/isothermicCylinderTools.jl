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
    return val1 * exp((im * w) * (th4prime(omega, tau)/th4(omega, tau)))
end

function rhombicAxisCalculation(w, omega, tau)
    val1 = im * th1prime(0, tau)/(2 * th2(omega, tau)) * th2(omega - im * w, tau)/th1(im * w, tau)
    return val1 * exp((im * w) * (th2prime(omega, tau)/th2(omega, tau)))
end



# VERY WIP
function numericallySolveRotation(v, wFunc, axisCalc)
    dwFunc = gradient(wFunc, v)[1]
    # ========== 2. Define right-hand side quaternion Q(v) ==========
    # Q(v) = sqrt(1 - w'(v)^2) * W1(w(v)) * k
    function Q_rhs(v)
        s = sqrt(1 - dwFunc(v)^2)
        wval = wFunc(v)
        W = axisCalc(wval)
        return Quaternion(0.0, 0.0, 0.0, s * W)  # Pure imaginary in k direction
    end

    # ========== 3. Define ODE: dΦ/dv = Q(v) * Φ(v) ==========
    function quat_ode!(dPhi, Phi_vec, v)
        # Convert Φ from vector to Quaternion
        Phi = Quaternion(Phi_vec...)
        Qv = Q_rhs(v)
        dQ = Qv * Phi
        # Convert back to vector form for the ODE solver
        dPhi[1] = dQ.s
        dPhi[2] = dQ.v[1]
        dPhi[3] = dQ.v[2]
        dPhi[4] = dQ.v[3]
    end

    # ========== 4. Solve the ODE numerically ==========
    # Initial condition: Φ(0) = 1 + 0i + 0j + 0k
    Phi0 = [1.0, 0.0, 0.0, 0.0]  # Represented as vector: [s, i, j, k]
    vspan = (0.0, 2pi)           # Solve from v = 0 to 2π

    prob = ODEProblem(quat_ode!, Phi0, vspan)
    return solve(prob, Tsit5())
end

using GLMakie
using Integrals
using Quaternions
using FiniteDifferences
using Zygote

"""
    integrateForm(A, B, x_0, y_0)

Integrates a pair of vector-valued differential forms `A` and `B` over a path from `(x_0, y_0)` to `(x, y)`.
Assumes the differential form is exact. Returns a function `f(x, y)` giving the immersion at `(x, y)`.

# Arguments
- `A`: Function `(x, y) -> (a1, a2, a3)` representing the first differential form.
- `B`: Function `(x, y) -> (b1, b2, b3)` representing the second differential form.
- `x_0`, `y_0`: Starting point for integration.

# Returns
- Function `(x, y) -> [X, Y, Z]` giving the integrated coordinates.
"""
function integrateForm(A, B, x_0, y_0)
    function f(x, y)
        gamma = (t) -> (x_0 + t*(x - x_0), y_0 + t*(y - y_0))

        dxdt = x-x_0
        dydt = y-y_0
        
        pullbackOmega1 = (t, p) -> A(gamma(t)...)[1] * dxdt + B(gamma(t)...)[1] * dydt        
        pullbackOmega2 = (t, p) -> A(gamma(t)...)[2] * dxdt + B(gamma(t)...)[2] * dydt
        pullbackOmega3 = (t, p) -> A(gamma(t)...)[3] * dxdt + B(gamma(t)...)[3] * dydt

        sol = []
        
        push!(sol, solve(IntegralProblem(pullbackOmega1, (0.0, 1.0)), QuadGKJL()).u)
        push!(sol, solve(IntegralProblem(pullbackOmega2, (0.0, 1.0)), QuadGKJL()).u)
        push!(sol, solve(IntegralProblem(pullbackOmega3, (0.0, 1.0)), QuadGKJL()).u)
        return sol
    end
    return f
end

"""
    numericDifferentials(f)

Computes the numerical partial derivatives of a vector-valued function `f(u, v)` with respect to `u` and `v`
using automatic differentiation.

# Arguments
- `f`: Function `(u, v) -> (x, y, z)`.

# Returns
- Tuple `(fu, fv)` where `fu(u, v)` and `fv(u, v)` are the partial derivatives with respect to `u` and `v`.
"""
function numericDifferentials(f)
    #fdm = central_fdm(n, d)  # 5-point stencil, 1st derivative
    
    #uF = (u, v) -> tuple([fdm(u -> f(u, v)[i], u) for i in 1:3]...)
    #vF = (u, v) -> tuple([fdm(v -> f(u, v)[i], v) for i in 1:3]...)

    uF = (u, v) -> tuple([something(Zygote.gradient(x -> f(x, v)[i], u)[1], 0.0) for i in 1:3]...)
    vF = (u, v) -> tuple([something(Zygote.gradient(x -> f(u, x)[i], v)[1], 0.0) for i in 1:3]...)
    return uF, vF
end

"""
    convertToQuaternion(f)

Converts a 3D vector-valued function `f(u, v)` to a quaternion-valued function with zero real part.

# Arguments
- `f`: Function `(u, v) -> (x, y, z)`.

# Returns
- Function `(u, v) -> Quaternion(0, x, y, z)`.
"""
function convertToQuaternion(f)
    return (u, v) -> Quaternions.Quaternion(
        0,
        f(u, v)[1], 
        f(u, v)[2], 
        f(u, v)[3]
    )
end

"""
    normalizeFunction(f)

Normalizes the output of a vector-valued function `f(u, v)` to unit length.

# Arguments
- `f`: Function `(u, v) -> (x, y, z)`.

# Returns
- Function `(u, v) -> (x', y', z')` where the vector has unit norm.
"""
function normalizeFunction(f)
    norm = (u, v) -> sqrt(f(u, v)[1]^2 + f(u, v)[2]^2 + f(u, v)[3]^2)
    return (u, v) -> (f(u, v)[1] / norm(u, v), f(u, v)[2] / norm(u,v), f(u, v)[3] / norm(u,v))
end

"""
    generateDualSurface(f; offsetX=0.0, offsetY=0.0)

Generates the dual surface of a parametric surface `f(u, v)` by integrating the normalized differentials.
Optionally offsets the base point.

# Arguments
- `f`: Parametric surface function `(u, v) -> (x, y, z)`.
- `offsetX`, `offsetY`: Optional offsets for the integration base point.

# Returns
- Function `(u, v) -> (x, y, z)` representing the dual surface.
"""
function generateDualSurface(f; offsetX = 0.0, offsetY = 0.0)
    fu, fv = numericDifferentials(f)

    outA = normalizeFunction(fu)

    outB = (u, v) -> (-1 * normalizeFunction(fv)(u,v)[1],
                      -1 * normalizeFunction(fv)(u,v)[2],
                      -1 * normalizeFunction(fv)(u,v)[3])
    
    return integrateForm(outA, outB, offsetX, offsetY)
end

"""
    generateDualSurface(fu, fv)

Generates normalized differential forms from given partial derivatives `fu` and `fv`.

# Arguments
- `fu`, `fv`: Functions representing partial derivatives with respect to `u` and `v`.

# Returns
- Tuple `(outA, outB)` of normalized differential forms.
"""
function generateDualSurface(fu, fv)
    outA = normalizeFunction(fu)

    outB = (u, v) -> (-1 * normalizeFunction(fv)(u,v)[1],
                      -1 * normalizeFunction(fv)(u,v)[2],
                      -1 * normalizeFunction(fv)(u,v)[3])
    
    return outA, outB
end

"""
    spinTransform(f, lambda)

Applies a quaternionic spin transformation to the surface `f(u, v)` using the quaternion-valued function `lambda(u, v)`.
Returns the integrated transformed surface.

# Arguments
- `f`: Parametric surface function `(u, v) -> (x, y, z)`.
- `lambda`: Function `(u, v) -> Quaternion`.

# Returns
- Function `(u, v) -> (x, y, z)` representing the transformed surface.
"""
function spinTransform(f, lambda)
    uF, vF = numericDifferentials(f)
    df1 = convertToQuaternion(uF)
    df2 = convertToQuaternion(vF)

    dfbar1 = (u, v) -> conj(lambda(u, v)) * df1(u, v) * lambda(u, v)
    dfbar2 = (u, v) -> conj(lambda(u, v)) * df2(u, v) * lambda(u, v)

    A = (u, v) -> (dfbar1(u, v).v1, dfbar1(u, v).v2, dfbar1(u, v).v3)
    B = (u, v) -> (dfbar2(u, v).v1, dfbar2(u, v).v2, dfbar2(u, v).v3)

    return integrateForm(A, B, 0.0, 0.0)
end

"""
    spinTransform(uF, vF, lambda)

Applies a quaternionic spin transformation to the differential forms `uF`, `vF` using `lambda(u, v)`.

# Arguments
- `uF`, `vF`: Functions representing partial derivatives.
- `lambda`: Function `(u, v) -> Quaternion`.

# Returns
- Tuple `(uF, vF)` of transformed differential forms.
"""
function spinTransform(uF, vF, lambda)
    df1 = convertToQuaternion(uF)
    df2 = convertToQuaternion(vF)

    dfbar1 = (u, v) -> conj(lambda(u, v)) * df1(u, v) * lambda(u, v)
    dfbar2 = (u, v) -> conj(lambda(u, v)) * df2(u, v) * lambda(u, v)

    uF = (u, v) -> (dfbar1(u, v).v1, dfbar1(u, v).v2, dfbar1(u, v).v3)
    vF = (u, v) -> (dfbar2(u, v).v1, dfbar2(u, v).v2, dfbar2(u, v).v3)

    return uF, vF
end
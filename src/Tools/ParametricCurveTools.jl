using GLMakie
using StaticArrays
using GeometryBasics
using Integrals
using FileIO
using MeshIO
using BSplineKit

"""
    plotParametricCurve(f, x; kwargs...)

Plot a parametric curve for parameter values `x`.

Arguments:
- `f::Function`: a function `f(t)` returning a length-3 tuple or vector (x,y,z).
- `x`: iterable of parameter values (e.g. `LinRange`, `range`, `Vector`).
- `kwargs...`: passed directly to Makie's `lines!`.

Returns the plot primitive returned by `lines!`.

Example:
    plotParametricCurve(t -> (sin(t), cos(t), t), LinRange(0, 2π, 200); color=:blue)
"""
function plotParametricCurve(f, x; kwargs...)
    lines!([[f(i)[k] for i in x] for k in 1:3]...; kwargs...)
end

"""
    plotCoordinateCurves(f, x, y, xCurveSamples, yCurveSamples; xColor = :red, yColor = :blue, kwargs...)

Plot coordinate (isoparametric) curves of a bivariate parametric surface `f(u,v)`.

Arguments:
- `f::Function`: function `f(u, v)` returning a length-3 tuple or vector (x,y,z).
- `x`, `y`: iterables of parameter values used to plot each curve (the sampling grid for the curve).
- `xCurveSamples`: iterable of `u` values for which `v`-curves (u fixed) are drawn.
- `yCurveSamples`: iterable of `v` values for which `u`-curves (v fixed) are drawn.
- `xColor`, `yColor`: colors for the two families of curves.
- `kwargs...`: forwarded to `plotParametricCurve` / Makie's `lines!`.

This function fixes one parameter at each sample and plots the parametric curve in the other parameter.
"""
function plotCoordinateCurves(f, x, y, xCurveSamples, yCurveSamples; xColor = :red, yColor = :blue, kwargs...   )
    for i in xCurveSamples
        plotParametricCurve((t) -> f(i, t), y; color = xColor, kwargs...)
    end
    for j in yCurveSamples
        plotParametricCurve((t) -> f(t, j), x; color = yColor, kwargs...)
    end
end


"""
    arcLengthFunction(curve) -> sFunc

Return a function `sFunc(t)` that computes the integral ∫_0^t norm(curve(u)) du using Integrals.jl.

Notes:
- The current integrand is `sqrt(curve(t)[1]^2 + curve(t)[2]^2)` (the Euclidean norm of the position vector),
  not the curve speed sqrt((dx/dt)^2 + (dy/dt)^2). If you need true arc length, replace the integrand with the
  magnitude of the derivative (e.g. compute derivatives with ForwardDiff).
- The returned `sFunc` evaluates the integral numerically on each call.
"""
function arcLengthFunction(curve)
    norm = (t, p) -> sqrt(curve(t)[1]^2 + curve(t)[2]^2)
    function outFunc(t)
        prob = IntegralProblem(norm, (0, t))
        sol = solve(prob, QuadGKJL())
        return sol.u
    end
    return outFunc
end

"""
    reparametrizeCurve(curve, tDomain, numSamples) -> (reparamCurve, total_s)

Build a reparameterized curve by sampling `curve` on `tDomain` and constructing a spline
interpolant from sampled arc-length values back to parameter `t`.

Arguments:
- `curve`: function `curve(t)` returning a length-2 or length-3 point.
- `tDomain`: tuple `(tmin, tmax, _)` where the first two entries give the parameter interval.
- `numSamples`: number of sample points to compute (used with `LinRange`).

Returns:
- `reparamCurve(s)`: a function taking `s` (arc-length-like parameter) and returning `curve(t)` evaluated
  at `t = U_of_S(s)` where `U_of_S` is the spline interpolant built from samples.
- `total_s`: the maximum sampled `s` value (useful as the slider upper bound).

Warnings:
- `sVal` must be monotonic for the inverse map `s -> t` to be well-defined.
- This function uses `arcLengthFunction`, which currently integrates position norm (see its docstring).
"""
function reparametrizeCurve(curve, tDomain, numSamples)
    sFunc = arcLengthFunction(curve)
    tVal = LinRange(tDomain[1], tDomain[2], numSamples)
    sVal = [sFunc(t) for t in tVal]
    U_of_S = interpolate(sVal, tVal, BSplineOrder(3))
    reparamCurve = (s) -> curve(U_of_S(s))
    return reparamCurve, maximum(sVal)
end
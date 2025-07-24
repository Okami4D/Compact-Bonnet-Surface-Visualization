using GLMakie
using StaticArrays
using GeometryBasics
using Integrals
using FileIO
using MeshIO

function plotParametricCurve(f, x; kwargs...)
    lines!([[f(i)[k] for i in x] for k in 1:3]...; kwargs...)
end

function plotCoordinateCurves(f, x, y, xCurveSamples, yCurveSamples; xColor = :red, yColor = :blue, kwargs...   )
    for i in xCurveSamples
        plotParametricCurve((t) -> f(i, t), y; color = xColor, kwargs...)
    end
    for j in yCurveSamples
        plotParametricCurve((t) -> f(t, j), x; color = yColor, kwargs...)
    end
end
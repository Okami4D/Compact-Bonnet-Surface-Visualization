using GLMakie
using StaticArrays
using GeometryBasics
using Integrals
using FileIO
using MeshIO

function plotParametricCurve(f, x; kwargs...)
    lines!([[f(i)[k] for i in x] for k in 1:3]...; kwargs...)
end

function plotCoordinateCurves(f, x, y, xCurveSamples, yCurveSamples; kwargs...)
    for i in xCurveSamples
        plotParametricCurve((t) -> f(i, t), y; color = :red)
    end
    for j in yCurveSamples
        plotParametricCurve((t) -> f(t, j), x; color = :blue)
    end
end
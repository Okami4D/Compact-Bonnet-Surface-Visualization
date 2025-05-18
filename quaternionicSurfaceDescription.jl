using GLMakie
using Quaternions
using Integrals
include("SurfaceViz.jl")

# Plot Setup
fig = Figure(size=(1920, 1080), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)

# Differential form df = A dx + B dy with f being the immersion of a sphere
A = (u, v) -> (-sin(v) * sin(u), sin(v) * cos(u), 0)
B = (u, v) -> (cos(v) * cos(u), cos(v) * sin(u), -sin(v))

#First version of integrating differential forms to obtain the corresponding immersion. Only works under the assumption
# that the differential form is exact.
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

f = integrateForm(A, B, 0.0, 0.0)
res = 50
phi = LinRange(0, -2π, res)
theta = LinRange(0, -π, res)

plotParametricSurface(f, phi, theta)
fig
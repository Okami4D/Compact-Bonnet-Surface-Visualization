using Integrals
using BSplineKit
using GLMakie

gamma = (t) -> (2 * sin(t), 2*cos(t))

function arcLengthFunction(curve)
    norm = (t, p) -> sqrt(curve(t)[1]^2 + curve(t)[2]^2)
    function outFunc(t)
        prob = IntegralProblem(norm, (0, t))
        sol = solve(prob, QuadGKJL())
        return sol.u
    end
    return outFunc
end

s = arcLengthFunction(gamma)


tDomain = (0, 2 * pi, 100)
tVal = LinRange(tDomain[1], tDomain[2], tDomain[3])

sVal = [s(t) for t in tVal]

U_of_S = interpolate(sVal, tVal, BSplineOrder(3))

gamma_Arc = (s) -> gamma(U_of_S(s))

fig = Figure(size = (1600, 1600))
ax1 = Axis(fig[1, 1])
xlims!(ax1, -3, 3)
ylims!(ax1, -3, 3)

slider = Slider(fig[2, 1], range = LinRange(0, maximum(sVal), 100), startvalue = 0)
on(slider.value) do t
    scatter!(ax1, [gamma_Arc(t)[1]], [gamma_Arc(t)[2]], color = :red)
    #scatter!(ax1, [gamma(t)[1]], [gamma(t)[2]], color = :blue)
end
fig

function reparametrizeCurve(curve, tDomain, numSamples)
    sFunc = arcLengthFunction(curve)
    tVal = LinRange(tDomain[1], tDomain[2], numSamples)
    sVal = [sFunc(t) for t in tVal]
    U_of_S = interpolate(sVal, tVal, BSplineOrder(3))
    reparamCurve = (s) -> curve(U_of_S(s))
    return reparamCurve, maximum(sVal)
end
using GLMakie
using Integrals
using Quaternions
using FiniteDifferences
using Zygote
using Alert
include("SurfaceViz.jl")

# Plot Setup
fig = Figure(size=(1500, 844), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)

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

function numericDifferentials(f)
    #fdm = central_fdm(n, d)  # 5-point stencil, 1st derivative
    
    #uF = (u, v) -> tuple([fdm(u -> f(u, v)[i], u) for i in 1:3]...)
    #vF = (u, v) -> tuple([fdm(v -> f(u, v)[i], v) for i in 1:3]...)

    uF = (u, v) -> tuple([something(Zygote.gradient(x -> f(x, v)[i], u)[1], 0.0) for i in 1:3]...)
    vF = (u, v) -> tuple([something(Zygote.gradient(x -> f(u, x)[i], v)[1], 0.0) for i in 1:3]...)
    return uF, vF
end

function convertToQuaternion(f)
    return (u, v) -> Quaternions.Quaternion(
        0,
        f(u, v)[1], 
        f(u, v)[2], 
        f(u, v)[3]
    )
end

function normalizeFunction(f)
    norm = (u, v) -> sqrt(f(u, v)[1]^2 + f(u, v)[2]^2 + f(u, v)[3]^2)
    return (u, v) -> (f(u, v)[1] / norm(u, v), f(u, v)[2] / norm(u,v), f(u, v)[3] / norm(u,v))
end

function generateDualSurface(f; offsetX = 0.0, offsetY = 0.0)
    fu, fv = numericDifferentials(f)

    outA = normalizeFunction(fu)

    outB = (u, v) -> (-1 * normalizeFunction(fv)(u,v)[1],
                      -1 * normalizeFunction(fv)(u,v)[2],
                      -1 * normalizeFunction(fv)(u,v)[3])
    
    return integrateForm(outA, outB, offsetX, offsetY)
end

function generateDualSurface(fu, fv)
    outA = normalizeFunction(fu)

    outB = (u, v) -> (-1 * normalizeFunction(fv)(u,v)[1],
                      -1 * normalizeFunction(fv)(u,v)[2],
                      -1 * normalizeFunction(fv)(u,v)[3])
    
    return outA, outB
end

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

function spinTransform(uF, vF, lambda)
    df1 = convertToQuaternion(uF)
    df2 = convertToQuaternion(vF)

    dfbar1 = (u, v) -> conj(lambda(u, v)) * df1(u, v) * lambda(u, v)
    dfbar2 = (u, v) -> conj(lambda(u, v)) * df2(u, v) * lambda(u, v)

    uF = (u, v) -> (dfbar1(u, v).v1, dfbar1(u, v).v2, dfbar1(u, v).v3)
    vF = (u, v) -> (dfbar2(u, v).v1, dfbar2(u, v).v2, dfbar2(u, v).v3)

    return uF, vF
end

f = parametricFuncPlane()

A = 1

fu, fv = numericDifferentials(f)
df1, df2 = generateDualSurface(fu, fv)


lamda1 = (u, v) -> Quaternions.Quaternion(0, A + f(u, v)[1], A + f(u, v)[2], A + f(u, v)[3])
lamda2 = (u, v) -> Quaternions.Quaternion(0, -A + f(u, v)[1], -A + f(u, v)[2], -A + f(u, v)[3])


gu1, gv1 = spinTransform(df1, df2, lamda1)
gu2, gv2 = spinTransform(df1, df2, lamda2)

g1 = integrateForm(gu1, gv1, 0.0, 0.0)
g2 = integrateForm(gu2, gv2, 0.0, 0.0)


res1 = 100

u1 = LinRange(-2, 2, res1)
v1 = LinRange(-2, 2, res1)

res2 = 10
u2 = LinRange(-2, 2, res2)
v2 = LinRange(-2, 2, res2)


plotParametricSurface(f, u1, v1)
plotParametricWireframe(g1, u2, v2)
plotParametricWireframe(g2, u2, v2)
alert("Julia script executed successfully!")
fig

#=
This Scene visualises the isothermic Torus with one family of curvatures lines given
by the paper by Bobenko and Hoffmann. We use "bobenkoCurves" and rotate this family to
generate a surface. The rotation is calculated by the function "numericallySolveRotation".

The specific example uses a tau-admissible reparametrization calculated numerically in the paper.
=#


using GLMakie
import Quaternions as Q

include("../../../Tools/ParametricSurfaceTools.jl")
include("../../../Tools/ParametricCurveTools.jl")
include("../../../Tools/isothermicCylinderTools.jl")

# Plot Setup
fig = Figure(
    size=(1200, 1200),
    scenekw=(
        lights=[DirectionalLight(RGBf(1, 1, 1),
            Vec3f(-1, 0, 0))],
    )
)

UVPlane = Axis(
    fig[1, 1:3],
    limits =(-2, 2, -2, 2),
    aspect = 1
)

LeftAx = Axis3(
    fig[2, 1],
    aspect=:data,
    perspectiveness=0.6,
    clip=false
)

MidAx = Axis3(
    fig[2, 2],
    aspect=:data,
    perspectiveness=0.6,
    clip=false
)

RightAx = Axis3(
    fig[2, 3],
    aspect=:data,
    perspectiveness=0.6,
    clip=false
)

hidedecorations!(LeftAx)
hidespines!(LeftAx)
hidedecorations!(MidAx)
hidespines!(MidAx)
hidedecorations!(RightAx)
hidespines!(RightAx)

fLeft = parametricFuncIsothermicHelicoid()

fMid = parametricFuncEnneper()

gamma = (t) -> (0.5 * sin(4 * t) + 1, t)
fRight, sMax = parametricFuncIsothermicSurfaceofRevolution(gamma, (0, 2*pi))


N = 100

umin, umax = -2, 2
vmin, vmax = -2, 2

u = LinRange(umin, umax, N)
v = LinRange(vmin, vmax, N)

plotParametricWireframe(fLeft, u, v, LeftAx; color = (:black, 0.05), transparency = true)
plotParametricWireframe(fMid, u, v, MidAx; color = (:black, 0.05), transparency = true)
plotParametricWireframe(fRight, u, v, RightAx; color = (:black, 0.05), transparency = true)

cursor_pos = Observable((0.0, 0.0))

on(events(UVPlane).mousebutton) do e
    if e.button == Mouse.left
        if is_mouseinside(UVPlane)
            x = mouseposition(UVPlane)[1]
            y = mouseposition(UVPlane)[2]
            cursor_pos[] = (x, y)
            notify!(cursor_pos)
        end
    end
end

deregister_interaction!(UVPlane,:rectanglezoom);

LineLen = 0.5


# Cross in 2D plot
xVals = lift(cursor_pos) do (u,v)
    points = [(u + i, v) for i in LinRange(-LineLen, LineLen, 20)]
end

yVals = lift(cursor_pos) do (u, v)
    points = [(u, v + i) for i in LinRange(-LineLen, LineLen, 20)]
end

lines!(UVPlane, xVals; color=:blue, linewidth=2)
lines!(UVPlane, yVals; color=:red, linewidth=2)

#function renderDownLines()
fLeftCrossX = lift(xVals) do X 
    points = [fLeft(p...) for p in X]
end
fLeftCrossY = lift(yVals) do Y 
    points = [fLeft(p...) for p in Y]
end

fMidCrossX = lift(xVals) do X 
    points = [fMid(p...) for p in X]
end
fMidCrossY = lift(yVals) do Y 
    points = [fMid(p...) for p in Y]
end

fRightCrossX = lift(xVals) do X 
    points = [fRight(p...) for p in X]
end
fRightCrossY = lift(yVals) do Y 
    points = [fRight(p...) for p in Y]
end

lines!(LeftAx, fLeftCrossX; color = :blue, linewidth = 2)
lines!(LeftAx, fLeftCrossY; color = :red, linewidth = 2)

lines!(MidAx, fMidCrossX; color = :blue, linewidth = 2)
lines!(MidAx, fMidCrossY; color = :red, linewidth = 2)

lines!(RightAx, fRightCrossX; color = :blue, linewidth = 2)
lines!(RightAx, fRightCrossY; color = :red, linewidth = 2)
#end




fig
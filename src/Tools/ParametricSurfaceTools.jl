using GLMakie
using StaticArrays
using GeometryBasics
using Integrals
using FileIO
using MeshIO

include("../Tools/ParametricCurveTools.jl")


"""
    plotParametricSurface(f, x, y; kwargs...)

Plots a parametric surface defined by the function `f(u, v)`, where `x` and `y` are ranges or vectors specifying the parameter domains for `u` and `v`. The function `f` should return a 3-element vector or tuple representing the (x, y, z) coordinates.

# Arguments
- `f`: Function of two variables `(u, v)` returning a 3-element vector or tuple.
- `x`: Range or vector for the first parameter.
- `y`: Range or vector for the second parameter.
- `kwargs...`: Additional keyword arguments passed to `surface`.

# Example
```julia
plotParametricSurface((u, v) -> [u, v, u^2 - v^2], -1:0.1:1, -1:0.1:1)
```
"""
function plotParametricSurface(f, x, y, axis::Union{Axis3, Nothing} = nothing; kwargs...)
    ax = isnothing(axis) ? current_axis() : axis
    surface!(ax, [[f(i,j)[k] for i in x, j in y] for k in 1:3]...; kwargs...)
end


"""
    plotParametricWireframe(f, x, y; kwargs...)

Plots a wireframe of a parametric surface defined by the function `f(u, v)`, where `x` and `y` are ranges or vectors specifying the parameter domains for `u` and `v`. The function `f` should return a 3-element vector or tuple representing the (x, y, z) coordinates.

# Arguments
- `f`: Function of two variables `(u, v)` returning a 3-element vector or tuple.
- `x`: Range or vector for the first parameter.
- `y`: Range or vector for the second parameter.
- `kwargs...`: Additional keyword arguments passed to `wireframe!`.

# Example
```julia
plotParametricWireframe((u, v) -> [u, v, sin(u)*cos(v)], 0:0.1:2π, 0:0.1:2π)
```
"""
function plotParametricWireframe(f, x, y, axis::Union{Axis3, Nothing} = nothing;  kwargs...)
    ax = isnothing(axis) ? current_axis() : axis
    wireframe!(ax, [[f(i,j)[k] for i in x, j in y] for k in 1:3]...; kwargs...)
    
end


"""
    generatePlanarBasedMesh(points, x, y, invertNormals=false)

Generates a mesh from a set of 3D points arranged in a grid defined by parameter ranges `x` and `y`. The mesh is constructed by connecting the points in a planar order, forming quadrilateral faces. Optionally, normals can be inverted.

# Arguments
- `points`: Vector of 3D points (e.g., `Point3f`) ordered according to the grid.
- `x`: Range or vector for the first parameter (number of columns).
- `y`: Range or vector for the second parameter (number of rows).
- `invertNormals`: (optional) If `true`, inverts the face normals.

# Returns
A `GeometryBasics.Mesh` object containing the mesh, normals, and color data.
"""
function generatePlanarBasedMesh(points, x, y, generateNormals = false)
    #= 
    Generate a list of faces, which consist of triangles whose vertecies 
    are specified by three indices in an array of type NgonFace.
    The indices are chosen in the natural order given by a plane.
    =#

    nx = length(x)
    ny = length(y)

    faces = [ 
        QuadFace(
            ((i +j*ny),
            (i+1 + j*ny),
            (i + 1 + (j+1)*ny),
            (i+(j+1)*ny))
        ) 
        for i in 1:(ny-1) for j in 0:(nx-2)
    ]

    #Calculate the Face Normals as additional information for lighting etc.
    if generateNormals 
        normalsCalc = face_normals(points, faces)
        return GeometryBasics.Mesh(
            points,
            faces,
            normal = normalsCalc,
        )
    else
        return GeometryBasics.Mesh(
            points,
            faces,
        )
    end
end


"""
    createParametricMesh(f, x, y; invertNormals=false)

Creates a mesh for a parametric surface defined by the function `f(u, v)`, where `x` and `y` are ranges or vectors specifying the parameter domains. The function `f` should return a 3-element vector or tuple representing the (x, y, z) coordinates for each `(u, v)`.

# Arguments
- `f`: Function of two variables `(u, v)` returning a 3-element vector or tuple.
- `x`: Range or vector for the first parameter.
- `y`: Range or vector for the second parameter.
- `invertNormals`: (optional) If `true`, inverts the face normals.

# Returns
A `GeometryBasics.Mesh` object representing the parametric surface.
"""
function createParametricMesh(f, x, y; invertNormals = false)
    points = Point3f[]
    for u in x, v in y
        p = Point3f(
            f(u, v)[1],
            f(u, v)[2],
            f(u, v)[3] 
        )
        push!(points, p)
    end
    return generatePlanarBasedMesh(points, x, y, invertNormals)
end


function rotateFunction(f, angle, axisVector)
    axis = Quaternion(0, axisVector...)
    rotQuat = Quaternion(cos(angle/2), 0, 0, 0) + sin(angle/2) * (1/sqrt(axis.s^2 + axis.v1^2 + axis.v2^2 + axis.v3^2)) * axis
    fRotatedQuaternion = (u, v) -> rotQuat * Quaternion(0, f(u, v)...) * (1/rotQuat)
    function output(u, v)
        fRotated = fRotatedQuaternion(u, v)
        return (fRotated.v1, fRotated.v2, fRotated.v3)
    end
    return output
end



#-------------------------------------------------------
#=
We define some basic embedding Functions which can be used in the above parametric plotting functions
=#
"""
    parametricFuncPlane()

Returns a parametric function for a plane in 3D.

# Returns
A function `(x, y) -> (x, y, 0)` representing the xy-plane.

# Example
```julia
plane = parametricFuncPlane()
x = -1:0.1:1
y = -1:0.1:1
# Use plane(x, y) in your plotting or mesh function
```
"""
function parametricFuncPlane()
    return (x, y) -> (x, y, 0)
end

"""
    parametricFuncSphere(r)

Returns a parametric function for a sphere of radius `r`.

# Arguments
- `r`: Radius of the sphere.

# Recommended Domains
For full coverage, use:
- `u ∈ [0, 2π]`
- `v ∈ [0, π]`

# Returns
A function `(u, v) -> (x, y, z)` representing the sphere.

# Example
```julia
sphere = parametricFuncSphere(1.0)
u = range(0, 2π, length=100)
v = range(0, π, length=100)
# Use sphere(u, v) in your plotting or mesh function
```
"""
function parametricFuncSphere(r)
    return (u, v) -> (
        r * sin(v) * cos(u),
        r * sin(v) * sin(u),
        r * cos(v)
    )
end

"""
    parametricFuncHelicoid(a, b)

Returns a parametric function for a helicoid surface.

# Arguments
- `a`: Radius scaling parameter.
- `b`: Height scaling parameter.

# Recommended Domains
For typical visualization, use:
- `u ∈ [-1, 1]`
- `v ∈ [0, 4π]`

# Returns
A function `(u, v) -> (x, y, z)` representing the helicoid.

# Example
```julia
helicoid = parametricFuncHelicoid(1.0, 0.2)
u = range(-1, 1, length=100)
v = range(0, 4π, length=100)
# Use helicoid(u, v) in your plotting or mesh function
```
"""
function parametricFuncHelicoid(a, b)
    return (u, v) -> (
        a * u * cos(v),
        a * u * sin(v),
        b * v
    )
end

"""
    parametricFuncTorus(r, R)

Returns a parametric function for a torus with minor radius `r` and major radius `R`.

# Arguments
- `r`: Minor (tube) radius.
- `R`: Major (center) radius.

# Recommended Domains
For full coverage, use:
- `u ∈ [0, 2π]`
- `v ∈ [0, 2π]`

# Returns
A function `(u, v) -> (x, y, z)` representing the torus.

# Example
```julia
torus = parametricFuncTorus(0.3, 1.0)
u = range(0, 2π, length=100)
v = range(0, 2π, length=100)
# Use torus(u, v) in your plotting or mesh function
```
"""
function parametricFuncTorus(r, R, domain = 2*pi, resolution = 100)
    return (u, v) -> (
        (r * cos(u) + R) * cos(v),
        (r * cos(u) + R) * sin(v),
        r * sin(u)
    )
end

"""
    parametricFuncEnneper()

Returns a parametric function for the Enneper minimal surface.

# Recommended Domains
For typical visualization, use:
- `u ∈ [-2, 2]`
- `v ∈ [-2, 2]`

# Returns
A function `(u, v) -> (x, y, z)` representing the Enneper surface.

# Example
```julia
enneper = parametricFuncEnneper()
u = range(-2, 2, length=100)
v = range(-2, 2, length=100)
# Use enneper(u, v) in your plotting or mesh function
```
"""
function parametricFuncEnneper()
    return (u, v) -> (
        u - (1/3) * u^3 + u * v^2,
        -v + (1/3) * v^3 - v * u^2,
        (u^2 - v^2)
    )
end


"""
    parametricFuncWente(theta; resolution=100)

Returns a parametric function for the Wente torus surface with parameter `theta`.

# Arguments
- `theta`: Parameter controlling the geometry of the Wente torus (in radians).

#Recommended Values
- `theta = 17.7324` (in degrees) is a common choice for the Wente torus.
- `theta = 12.7898` (in degrees) is another option.
- `theta = 21.4807` (in degrees) is also used.
- `theta = 9.9285` (in degrees) is a less common choice.
- `theta = 22.8449` (in degrees) is a less common choice.

# Recommended Domains
For typical visualization, use:
- `u ∈ [-π/2, π/4]`
- `v ∈ [-π, 2π/3]`

# Returns
A function `(u, v) -> (x, y, z)` representing the Wente torus surface, suitable for use with parametric plotting or mesh generation routines.

# Example
```julia
wente = parametricFuncWente(1.2)
u = range(-π/2, π/4, length=100)
v = range(-π, 2π/3, length=100)
# Use wente(u, v) in your plotting or mesh function
```
"""
function parametricFuncWente(thetaDeg = 17.7324)
    theta = thetaDeg * pi / 180
    H = 1/2
    thetaBar = (65.35 * pi) / 180

    k = sin(theta)
    kbar = sin(thetaBar)
    gamma = sqrt(tan(theta))
    gammabar = sqrt(tan(thetaBar))

    alpha = sqrt((4 * H * sin(2 * thetaBar))/(sin(2 * (theta + thetaBar))))
    alphabar = sqrt((4 * H * sin(2 * theta))/(sin(2 * (theta + thetaBar))))

    b = - (4 * H * sin(2 * theta) * cos( 2 * thetaBar))/ (sin(2 * (theta + thetaBar)))
    p = (4 * H * sin(2 * theta) * sin(2 * thetaBar))/(sin(2 * (theta + thetaBar)))

    Gamma = gamma * gammabar

    Z(u, v) = (sqrt(2) * ((( alphabar^2 - b)*(gamma * cos(u))^2 + p) * gammabar *  cos(v) - (p * (gamma * cos(u))^2 + alphabar^2 + b) * gamma * cos(u))) / (( sqrt(H) * alphabar^2) * ((1 - Gamma * cos(u)*cos(v)) * sqrt(p - 2 * b * (gamma * cos(u))^2 - p * (gamma * cos(u))^4)))
    
    w(u) = (1/alpha) * 2 * sqrt(H) * solve(IntegralProblem((t, p) -> (1 + (Gamma * cos(t))^2)/((1 - (Gamma * cos(t))^2) * sqrt(1 - (k * sin(t))^2)), (0.0, u)), QuadGKJL()).u

    j(u) = atan((tan(u) * alpha) * sqrt(1 - k^2 * sin(u)^2 )/(2 * sqrt(H))) * ceil(trunc(2*u / pi)/2) * pi

    I(v) = solve(IntegralProblem((t, p) -> (1 - 2 * (kbar * sin(t))^2)/(sqrt(1 - (kbar * sin(t))^2)), (0.0, v)), QuadGKJL()).u
    
    x(u, v) = Z(u, v) * cos(w(u) - j(u)) + cos(w(u)/2 * H)
    y(u, v) = Z(u, v) * sin(w(u) - j(u)) + sin( w(u)/2 * H)
    z(u, v) = 1/(alphabar * sqrt(H)) * ((2 * Gamma * cos(u) * sin(v) * sqrt(1 - (kbar * sin(v))^2 ))/(1 - Gamma * cos(u) * cos(v)) .+ (1/gammabar) * I(v))
    return (u, v) -> (x(u, v), y(u, v), z(u, v))
end

"""
    parametricFuncKleinBottle(a=2.0, b=1.0)

Returns a parametric function for the Klein bottle surface.

# Arguments
- `a`: (optional) Main radius of the tube (default: 2.0).
- `b`: (optional) Radius of the tube (default: 1.0).

# Recommended Domains
For typical visualization, use:
- `u ∈ [0, 2π]`
- `v ∈ [0, 2π]`

# Returns
A function `(u, v) -> (x, y, z)` representing the Klein bottle.

# Example
```julia
klein = parametricFuncKleinBottle()
u = range(0, 2π, length=100)
v = range(0, 2π, length=100)
# Use klein(u, v) in your plotting or mesh function
```
"""
function parametricFuncKleinBottle(a=2.0, b=1.0)
    return (u, v) -> begin
        x = (a + b * cos(u/2) * sin(v) - b * sin(u/2) * sin(2v)) * cos(u)
        y = (a + b * cos(u/2) * sin(v) - b * sin(u/2) * sin(2v)) * sin(u)
        z = b * sin(u/2) * sin(v) + b * cos(u/2) * sin(2v)
        (x, y, z)
    end
end

"""
    parametricFuncCylinder(r, h)
Returns a parametric function for a cylinder of radius `r` and height `h`.
    This is isothermic for r = 1.
# Arguments
- `r`: Radius of the cylinder.
# Recommended Domains
For typical visualization, use:
- `u ∈ [0, 2π]`
- `v ∈ [0, height]`
# Returns
A function `(u, v) -> (x, y, z)` representing the cylinder.
# Example
```julia
cylinder = parametricFuncCylinder(1.0)
u = range(0, 2π, length=100)
v = range(0, 5.0, length=100)
# Use cylinder(u, v) in your plotting or mesh function
``` 
"""
function parametricFuncCylinder(r)
    return (u, v) -> (
        r * cos(u),
        r * sin(u),
        v
    )
end

function parametricFuncIsothermicCatenoid()

    out = (theta, phi) -> begin
        t = atanh(tan(phi/2)) * 2
        x = cosh(t) * cos(theta)
        y = cosh(t) * sin(theta)
        z = t
        return (x, y, z)
    end
    return out
end

function parametricFuncIsothermicCone(k)
    
    out = (theta, psi) -> begin
        s = exp(psi/sqrt(1 + k^2))
        x = s * cos(theta)
        y = s * sin(theta)
        z = k * s
        return (x, y, z)
    end
    return out
end

function parametricFuncIsothermicSurfaceofRevolution(gamma, domain; nSamples = 100)
    reparam = reparametrizeCurve(gamma, domain, nSamples)
    gamma_Arc = reparam[1]
    out = (theta, psi) -> begin
        r, z = gamma_Arc(psi)
        x = r * cos(theta)
        y = r * sin(theta)
        return (x, y, z)
    end
    return out, reparam[2]
end


function parametricFuncModulatedTorus(modFunc, R)
    return (phi, theta) -> begin
        mod = modFunc(phi, theta)
        x = ((R + mod * cos(theta)) * cos(phi))
        y = ((R + mod * cos(theta)) * sin(phi))
        z = mod * sin(theta)
        return (x, y, z)
    end
end

function parametricFuncIsothermicHelicoid()
    return (v, u) -> (sinh(u) * cos(v), sinh(u) * sin(v), v)
end
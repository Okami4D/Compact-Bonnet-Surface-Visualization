using GLMakie
using StaticArrays
using GeometryBasics
using Integrals
using FileIO
using MeshIO

#-------------------------------------------------------
#=
This method generates a mesh from a set of points
and a resolution. It creates a mesh of triangles assuming
that all points are given in the order of a plane.
=#
function generatePlanarBasedMesh(points, resolution = 100, invertNormals = false)
    #= 
    Generate a list of faces, which consist of triangles whose vertecies 
    are specified by three indices in an array of type NgonFace.
    The indices are chosen in the natural order given by a plane.
    =#
    faceList = NgonFace{3, Int64}[]

    for i in 1:resolution-1, j in 1:resolution-1
        idx = (j - 1) * resolution + i

        p1 = idx
        p2 = idx + resolution
        p3 = idx + 1
        p4 = idx + resolution + 1

        if invertNormals
            push!(faceList, NgonFace{3, Int64}([p1, p3, p2]))
            push!(faceList, NgonFace{3, Int64}([p2, p3, p4]))
        else
            push!(faceList, NgonFace{3, Int64}([p1, p2, p3]))
            push!(faceList, NgonFace{3, Int64}([p3, p2, p4]))
        end



    end

    #Calculate the Face Normals as additional information for lighting etc.
    normalsCalc = face_normals(points, faceList)

    #Output the entirty of information so far as a mesh object from GeometryBasics to bundle it all together
    outputMesh = GeometryBasics.Mesh(
        points,
        faceList,
        normal = normalsCalc,
        #We hardcode the color for now, but this could be changed to a gradient or something else
        color = [v[3] for v in points]
    )

    return outputMesh
end

#=
This method generates a mesh given a set of parametric functions
    - xfunc: function for x coordinate
    - yfunc: function for y coordinate
    - zfunc: function for z coordinate
    - u: list of u values
    - v: list of v values

Keep in mind that u and v must be of the same length
    - resolution: number of sampling points in each direction
=#
function createParametricSurface(xfunc, yfunc, zfunc, u, v, invertNormals = false)
    if length(u) != length(v)
        throw(ArgumentError("u and v must have the same length for a square grid"))
    end

    points = Point3f[]
    for uu in u, vv in v
        # Compute vertex indices (row-major order)
        p1 = Point3f(
            xfunc(uu, vv), 
            yfunc(uu, vv), 
            zfunc(uu, vv)
        )
        push!(points, p1)
    end
    return generatePlanarBasedMesh(points, length(u), invertNormals)
end


#-------------------------------------------------------
#=
We define some functions to create some basic shapes for testing
and debugging purposes.

List of Objects created:
    - Plane
    - Sphere
    - Helicoid
    - Torus
=#


#=
Generate a plane with dimensions 2a x 2b
    - a: half the width of the plane
    - b: half the height of the plane
    - resolution: number of sampling points in each direction
=#
function createPlane(a, b, resolution=100)

    u = range(-a, a, length=resolution)
    v = range(-b, b, length=resolution)

    x(u, v) = u
    y(u, v) = v
    z(u, v) = 0

    return createParametricSurface(x, y, z, u, v)
end


#=
Generate a sphere with radius r
    - r: radius of the sphere
    - resolution: number of sampling points in each direction
=#
function createSphere(r, resolution = 100)

    u = range(0, 2π, length=resolution)
    v = range(0, π, length=resolution)

    x(u, v) = r*sin(v)*cos(u)
    y(u, v) = r*sin(v)*sin(u)
    z(u, v) = r*cos(v)
    return createParametricSurface(x, y, z, u, v, true)
end

#=
Generate a helicoid with parameters a and b
    - a: radius of the helicoid
    - b: height of the helicoid
    - resolution: number of sampling points in each direction
=#
function createHelicoid(a, b, domain = 20, resolution = 100)
    u = range(0, domain, length=resolution)
    v = range(0, domain, length=resolution)

    x(u, v) = a*u*cos(v)
    y(u, v) = a*u*sin(v)
    z(u, v) = b*v
    return createParametricSurface(x, y, z, u, v)
end

#=
Generate a Torus with parameters r and R
    - r: Small Circle raidus - thickness
    - r: big Radius
    - domain: range of u and v
    - resolution: number of sampling points in each direction
=#
function createTorus(r, R, domain = 2*pi, resolution = 100)
    u = range(0, domain, length=resolution)
    v = range(0, domain, length=resolution)

    x(u, v) = (r*sin(u) + R) * cos(v)
    y(u, v) = (r*sin(u) + R) * sin(v)
    z(u, v) = r*cos(u)
    return createParametricSurface(x, y, z, u, v, false)
end

#=
Generate an Enneper surface
    - domain: range of u and v
    - resolution: number of sampling points in each direction
=#
function createEnneper(domain = 5, resolution = 100)
    u = range(-domain, domain, length=resolution)
    v = range(-domain, domain, length=resolution)

    x(u, v) = u - (1/3) * u^3 + u * v^2
    y(u, v) = -v + (1/3) * v^3 - v * u^2
    z(u, v) = (u^2 - v^2)
    return createParametricSurface(x, y, z, u, v, false)
end

# TODO: Add a function to create a Wente surface
function createWente(theta, domainU = 2.761094, domainV = 3.58655, resolution = 40)
    u_ = range(-pi /2, domainU, length=resolution)
    v_ = range(-pi, domainV, length=resolution)

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
    return createParametricSurface(x, y, z, u_, v_, false)
end

#-------------------------------------------------------
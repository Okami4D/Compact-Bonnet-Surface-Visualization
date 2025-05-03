using GLMakie
using StaticArrays
using GeometryBasics
using FileIO
using MeshIO

#-------------------------------------------------------
#=
This method generates a mesh from a set of points
and a resolution. It creates a mesh of triangles assuming
that all points are given in the order of a plane.
=#
function generatePlanarBasedMesh(points, resolution = 100)
    #= 
    Generate a list of faces, which consist of triangles whose vertecies 
    are specified by three indices in an array of type NgonFace.
    The indices are chosen in the natural order given by a plane.
    =#
    faceList = NgonFace{3, Int64}[]

    for i in 1:resolution-1, j in 1:resolution-1
        idx = (j - 1) * resolution + i
        push!(faceList, NgonFace{3, Int64}([idx + 1, idx + resolution, idx + resolution + 1]))
        push!(faceList, NgonFace{3, Int64}([idx, idx + resolution, idx + 1]))
        
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
function createParametricSurface(xfunc, yfunc, zfunc, u, v)
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
    return generatePlanarBasedMesh(points, length(u))
end


#-------------------------------------------------------
#=
We define some functions to create some basic shapes for testing
and debugging purposes.

List of Objects created:
    - Plane
    - Sphere
    - Helicoid
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
    return createParametricSurface(x, y, z, u, v)
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

#-------------------------------------------------------
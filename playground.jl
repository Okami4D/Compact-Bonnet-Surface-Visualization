include("SurfaceViz.jl")

# Plot Setup
fig = Figure(size=(1920, 1080), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.8, 
    clip=false
)
hidedecorations!(ax)
hidespines!(ax)






#mesh!(ax, createPlane(2, 2))
#mesh!(ax, createSphere(2))

hel = createHelicoid(2, 10, 60, 400)
mesh!(ax, hel, specular = 0.4, diffuse = 0.7)

#save("helicoid.obj", hel)

fig
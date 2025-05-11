using GLMakie
include("SurfaceViz.jl")

# Plot Setup
fig = Figure(size=(1920, 1080), scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :equal, 
    perspectiveness = 0.8, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)

activeMesh = []

# Generate Mesh 
activeMesh = createEnneper(5, 40)
mesh!(ax, activeMesh, specular = 0.4, diffuse = 0.7)


# Setup the Buttons
fig[2, 1] = buttongrid = GridLayout(tellwidth = false)
buttons = buttongrid[1, 1:4] = [
    Button(fig, label = "Surface", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Wireframe", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Save PNG", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Save OBJ", width = 150, height = 70, fontsize = 30)
]

on(buttons[1].clicks) do b
    empty!(ax)
    mesh!(ax, activeMesh, specular = 0.4, diffuse = 0.7)
    fig
end
on(buttons[2].clicks) do b
    empty!(ax)
    wireframe!(ax, activeMesh, linewidth = 0.5)
    fig
end
on(buttons[3].clicks) do b
    save("output.png", fig, size = (1920, 1080))
    fig
end
on(buttons[4].clicks) do b
    save("output.obj", hel)
    fig
end

fig
#=
#Rotating Camera around center
fps = 60
nframes = 1200

for i in 1:nframes
    ang[] = ang[] + 0.01
    sleep(1 / fps)
end
#save("helicoid.obj", hel)
=#
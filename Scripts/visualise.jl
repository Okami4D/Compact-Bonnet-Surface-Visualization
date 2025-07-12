using GLMakie
include("Tools/ParametricSurfaceTools.jl")

# Plot Setup
fig = Figure(
    size=(1200, 800), 
    scenekw = (lights = [
        DirectionalLight(RGBf(1, 1, 1), 
        Vec3f(-1, 0, 0))
    ],
))

ang = Observable(1.01 * pi)

ax = Axis3(
    fig[1, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)


# Generate Mesh 
resolution = 100
x = LinRange(0, 2 * pi, resolution)
y = LinRange(0, 2 * pi, resolution)
#activeParametrization = parametricFuncWente(9.9285)
activeParametrization = parametricFuncKleinBottle()
#activeParametrization = parametricFuncEnneper()
activeMesh = createParametricMesh(activeParametrization, x, y)

plotParametricSurface(activeParametrization, x, y;
        specular = 0.4,
        diffuse = 0.7
) 


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
    plotParametricSurface(activeParametrization, x, y; 
        specular = 0.4, 
        diffuse = 0.7
    )    
end
on(buttons[2].clicks) do b
    empty!(ax)
    plotParametricWireframe(activeParametrization, x, y)
end
on(buttons[3].clicks) do b
    save("output.png", fig, update = false)
end
on(buttons[4].clicks) do b
    save("output.obj", activeMesh)
end

fig
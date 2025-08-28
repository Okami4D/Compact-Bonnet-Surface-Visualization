using GLMakie
include("../../Tools/ParametricSurfaceTools.jl")

resolution = 120


# Plot Setup
fig = Figure(
    size=(1200, 800), 
    scenekw = (lights = [
        DirectionalLight(RGBf(1, 1, 1), 
        Vec3f(-1, 0, 0))
    ],
))

ax = Axis3(
    fig[2, 1], 
    aspect = :data, 
    perspectiveness = 0.6, 
    clip=false
)

hidedecorations!(ax)
hidespines!(ax)




topGrid = GridLayout(fig[1, 1], rows = 1, columns = 2, tellwidth = false)
menu = Menu(topGrid[1, 1], 
    options = ["Plane", "Sphere", "Helicoid", "Torus", "Enneper", "Wente", "Klein Bottle"],
    default = "Plane"
)


# Setup the Buttons
fig[3, 1] = buttongrid = GridLayout(tellwidth = false)
buttons = buttongrid[1, 1:4] = [
    Button(fig, label = "Surface", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Wireframe", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Save PNG", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Save OBJ", width = 150, height = 70, fontsize = 30)
]
# Set up the event listeners for the buttons
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


on(menu.selection) do selection
    if selection == "Plane"
        render!(parametricFuncPlane(), -2, 2, -2, 2, resolution)
    elseif selection == "Sphere"
        render!(parametricFuncSphere(1), 0, pi, 0, 2 * pi, resolution)
    elseif selection == "Helicoid"
        render!(parametricFuncHelicoid(1, 1), -1, 1, 0, 4 * pi, resolution)
    elseif selection == "Torus"
        render!(parametricFuncTorus(1, 2), 0, 2 * pi, 0, 2 * pi, resolution)
    elseif selection == "Enneper"
        render!(parametricFuncEnneper(), -2, 2, -2, 2, resolution)
    elseif selection == "Wente"
        render!(parametricFuncWente(17.7324), -pi/2, pi/4, -pi, 2 * pi/3, resolution)
    elseif selection == "Klein Bottle"
        render!(parametricFuncKleinBottle(), 0, 2 * pi, 0, 2 * pi, resolution)
    else
        error("Unknown selection: $selection")
    end
end




function render!(parametrization, xmin, xmax, ymin, ymax, resolution)
    empty!(ax)
    x = LinRange(xmin, xmax, resolution)
    y = LinRange(ymin, ymax, resolution)
    
    plotParametricSurface(parametrization, x, y;
            specular = 0.4,
            diffuse = 0.7
    )
    fig
end

render!(parametricFuncPlane(), -2, 2, -2, 2, resolution)
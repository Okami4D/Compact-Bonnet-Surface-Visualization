using GLMakie
include("../../Tools/ParametricSurfaceTools.jl")

global renderParam = parametricFuncPlane()
global renderRes = 50
global renderX = LinRange(-2, 2, renderRes)
global renderY = LinRange(-2, 2, renderRes)
global wireframeRendering = false



# Plot Setup
fig = Figure(
    size=(1500, 1000), 
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
ObjectMenu = Menu(topGrid[1, 1], 
    options = ["Plane", "Sphere", "Helicoid", "Torus", "Enneper", "Wente", "Klein Bottle", "Isothermic Catenoid", "Isothermic Cone"],
    default = "Plane"
)

renderingMenu = Menu(topGrid[1, 2], 
    options = ["Surface", "Wireframe"],
    default = "Surface"
)


# Setup the Buttons
fig[3, 1] = buttongrid = GridLayout(tellwidth = false)
buttons = buttongrid[1, 1:2] = [
    Button(fig, label = "Save PNG", width = 150, height = 70, fontsize = 30),
    Button(fig, label = "Save OBJ", width = 150, height = 70, fontsize = 30)
]

resolutionSlider = Slider(fig[4, 1], range = 10:10:200, startvalue = 50)

on(resolutionSlider.value) do v
    global renderRes = Int(v)
    global renderX = LinRange(minimum(renderX), maximum(renderX), renderRes)
    global renderY = LinRange(minimum(renderY), maximum(renderY), renderRes)
    render!()
end

on(buttons[1].clicks) do b
    save("output.png", fig, update = false)
end
on(buttons[2].clicks) do b
    save("output.obj", createParametricMesh(renderParam, renderX, renderY))
end



on(renderingMenu.selection) do selection
    if selection == "Surface"
        global wireframeRendering = false
        render!()
    elseif selection == "Wireframe"
        global wireframeRendering = true
        render!()
    end
end



on(ObjectMenu.selection) do selection
    if selection == "Plane"
        render!(parametricFuncPlane(), -2, 2, -2, 2)
    elseif selection == "Sphere"
        render!(parametricFuncSphere(1), 0, pi, 0, 2 * pi)
    elseif selection == "Helicoid"
        render!(parametricFuncHelicoid(1, 1), -1, 1, 0, 4 * pi)
    elseif selection == "Torus"
        render!(parametricFuncTorus(1, 2), 0, 2 * pi, 0, 2 * pi)
    elseif selection == "Enneper"
        render!(parametricFuncEnneper(), -2, 2, -2, 2)
    elseif selection == "Wente"
        render!(parametricFuncWente(17.7324), -pi/2, pi/4, -pi, 2 * pi/3)
    elseif selection == "Klein Bottle"
        render!(parametricFuncKleinBottle(), 0, 2 * pi, 0, 2 * pi)
    elseif selection == "Isothermic Catenoid"
        render!(parametricFuncIsothermicCatenoid(), 0, 2 * pi, -pi/2 + 0.3,  pi/2 - 0.3)
    elseif selection == "Isothermic Cone"
        render!(parametricFuncIsothermicCone(1), 0, 2 * pi, -1, 1)
    else
        error("Unknown selection: $selection")
    end
end




function render!()
    empty!(ax)
    if wireframeRendering
        plotParametricWireframe(renderParam, renderX, renderY;  color = (:black, 0.1), transparency = true)
    else
        plotParametricSurface(renderParam, renderX, renderY; 
            specular = 0.4, 
            diffuse = 0.7
        )
    end
    fig
end

function render!(parametrization, xmin, xmax, ymin, ymax)
    global renderParam, renderRes, renderX, renderY
    renderParam = parametrization
    renderX = LinRange(xmin, xmax, renderRes)
    renderY = LinRange(ymin, ymax, renderRes)
    render!()
end

render!()
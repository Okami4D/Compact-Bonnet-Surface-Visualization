This Project aims to visualise constructions given in a paper by Alexander Bobenko, Tim Hoffmann and Andrew O. Sageman-Furnas called [“Compact Bonnet Pairs: isometric tori with the same curvatures”](https://arxiv.org/abs/2110.06335). It was created under the supervision of Prof. Tim Hoffmann at the TUM School of Computation, Information and Technology. 

The central aim of the paper was to use constructions from the Paper [“Isothermic Tori with one family of planar curvature lines and area contrained hyperbolic elastica”](https://arxiv.org/abs/2312.14956) to explicitly construct pairs of immersed real analytic tori in three dimensional Euclidean space. This project applies the techniques introduced and implements them into interactive scenes using GLMakie. 

We divide up the scripts used into “Tools” and “Scripts”. The “Tools” build the programmatic framework for rendering and implementing the general theory given in the papers. Here is a quick list of functions covered in the different scripts:

1. “Parametric Curve Tools” - Functions (mostly Wrappers for GLMakie) allowing us to render 2D-Curves
2. “Parametric Surface Tools” - Functions (partly Wrappers for GLMakie and partly custom maps) allowing us to render general parametric surfaces using a map ``F: \mathbb{R}^2 \to \mathbb{R}^3``. It also contains a list of pre-defined parametrizations to test with.
3. “Quaternionic Geometry Toolkit” - Formalization of the quaternionic description of immersed surfaces. Explicit generation of christoffel dual surfaces, spin transforms and numerical integration of a given differential form.
4. “isothermic Cylinder Tools” - Implementation of the “Isothermic Tori” paper with numerical computation to generate the immersion of an isothermic cylinder given ``w, \tau, \omega``
5. “jacobi Theta Extension” - An extension of the jacobi theta functions used in the project extending the functionality slightly from the EllipticFunctions.jl package

## Visualization Baseline

The framework constructed is aimed to visualize any parametric function
```math
f: U \to \mathbb{R}^{3}
```
where ``U = [x_{min}, x_{max}] \times [y_{min}, y_{max}] \subset \mathbb{R}^{3}``.

This base uses GLMakie for rendering and simply creates wrappers around the “surface!” and “wireframe!” functions to allow inputting literal maps $f$ and $x$ and $y$ domain values generating $U$.

There is also a simple construction of a GeometryBasics Mesh given the “createParametricMesh” function which allows for outputting the parametrization of our function into an “.obj” for futher visualization in programms like Blender (This functionality can be seen in the Interactive Scene “visualiseParametricSurfaces.jl”)

Eventhough the Project centers around generating parametric functions for very complicated constructions (namely the isothermic Torus and its Bonnet Surfaces) we have added some simple pre-built parametric functions that can all be seen in action in the “visualiseParametricSurfaces.jl” scene.

One note is that the parametrization given for the Wente Torus, a closely related surface to the topic of isothermic Tori and Bonnet Pairs is taken from Rolf Walter’s “Explicit Examples to the H-Problem of Heinz Hopf”.

A very similar wrapper structure is only lightly implemented for curves in “ParametricCurveTools.jl”. This is only used in one Stock interactive scene (namely in “Bobenko Curves”) but can still be just as useful as the wrappers for surfaces. The reason for creating these is just to unify the language used when creating a scene to actually work with parametric function $f$ and avoid explicitly defining points of a surface.

## Quaternionic Geometry

In the main paper a known framework for working with three dimensional surfaces as quaternions is introduced. It allows some smart and fast calculations especially if we want to express the Bonnet Problem explicitly. We implement aspects of this model by using “Quaternions.jl” in the file “QuaternionicGeometryToolkit.jl”.

The main aim of this tool collection is to take a given parametrization and make it possible to apply either the spin-transform or calculate the christoffel dual of the surface. These operations are given by considering a parametrization

```math
f: \mathbb{R}^{2} \to \mathbb{R}^{3}, \quad \begin{pmatrix}x \\ y\end{pmatrix} \mapsto \begin{pmatrix}f_{1}(x, y ) \\ f_{2} (x, y) \\ f_{3}(x, y)\end{pmatrix}
```

as a map

```math
\tilde{f}: \mathbb{R}^{2} \to \mathbb{H}, \quad \begin{pmatrix}x \\ y\end{pmatrix} \mapsto f_{1}(x, y) \mathbb{i} + f_{2}(x, y) \mathbb{j} + f_{3}(x, y) \mathbb{k}
```

this is achieved programmatically by the function “convertToQuaternion”. We can then consider the differential of this map 
```math
d \tilde{f} = \del_{x} f dx + \del_{y} f dy
```

Then the operations we are concerend with are given for the
- **Christoffel Dual**: We call a parametrization ``h: \mathbb{R}^{2} \to \mathbb{R}^{3}`` the christoffel dual of $f$ if 
```math
dh = e^{-2h} (\del_{u}f du - \del_{v}f dv) = \frac{\del_{u}f}{\norm{\del_{u}f}}du - \frac{\del_{v}f}{\norm{\del_{v}f}}dv
```


- **Spin Transform**: Given a quaternion ``\lambda: \mathbb{R}^{2} \to \mathbb{H}`` we define the spin transform of ``f`` by ``\lambda`` to be the parametrization $h$ satisfying: 

```math
dh = \overline{\lambda} df \lambda
```

Because all of these definitions only define the differential of our parametrizations the toolkit also provides a compatible way to numerically integrate the parametrization. This is not very performant since in the aim to keep the outputs of all operations functions that act as parametrizations, we require that at every evaluation we recalculate the numerical integration and differentiation. Meaning that for every $(x, y)$ we have the sequence
```math
(x, y) \rightarrow \text{numerically differentiate} \rightarrow \text{evaluate operation} \rightarrow \text{numerically integrate} \rightarrow \begin{pmatrix}X \\ Y \\Z\end{pmatrix}
```
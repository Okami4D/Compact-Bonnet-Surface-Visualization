This Project was built under the supervision of Prof. Dr. Tim Hoffmann at the Technical University of Munich. It aims to visualize the paper "Compact Bonnet Pairs: Isometric Tori with the same Curvatures" and the applied techniques in interactive GLMakie scenes.

Math Test
``\det(\lambda I_{k} - \mathcal{J}_{k}) = \begin{cases}(\lambda^{2}-1)^{\frac{k}{2}} &\quad k \text{ even} \\(\lambda - 1)(\lambda^{2}-1)^{\frac{k-1}{2}}  &\quad k \text{ odd}\end{cases}``

This Project aims to visualise constructions given in a paper by Alexander Bobenko, Tim Hoffmann and Andrew O. Sageman-Furnas called “Compact Bonnet Pairs: isometric tori with the same curvatures”. It was created under the supervision of Prof. Tim Hoffmann at the TUM School of Computation, Information and Technology. 

The central aim of the paper was to use constructions from the Paper “Isothermic Tori with one family of planar curvature lines and area contrained hyperbolic elastica” to explicitly construct pairs of immersed real analytic tori in three dimensional Euclidean space. This project applies the techniques introduced and implements them into interactive scenes using GLMakie. 

We divide up the scripts used into “Tools” and “Scripts”. The “Tools” build the programmatic framework for rendering and implementing the general theory given in the papers. Here is a quick list of functions covered in the different scripts:

1. “Parametric Curve Tools” - Functions (mostly Wrappers for GLMakie) allowing us to render 2D-Curves
2. “Parametric Surface Tools” - Functions (partly Wrappers for GLMakie and partly custom maps) allowing us to render general parametric surfaces using a map ``F: \mathbb{R}^2 \to \mathbb{R}^3``. It also contains a list of pre-defined parametrizations to test with.
3. “Quaternionic Geometry Toolkit” - Formalization of the quaternionic description of immersed surfaces. Explicit generation of christoffel dual surfaces, spin transforms and numerical integration of a given differential form.
4. “isothermic Cylinder Tools” - Implementation of the “Isothermic Tori” paper with numerical computation to generate the immersion of an isothermic cylinder given ``w, \tau, \omega``
5. “jacobi Theta Extension” - An extension of the jacobi theta functions used in the project extending the functionality slightly from the EllipticFunctions.jl package
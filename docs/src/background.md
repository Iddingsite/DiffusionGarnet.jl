# [Background](@id background)

### Theory

Garnet is a mineral commonly used in metamorphic petrology to better understand geological processes, as it occurs in a variety of rock types. This mineral often exhibits a wide range of compositional zoning, which can be interpreted as recording ranges of pressure and temperature (PT) conditions. Modelling diffusion processes can help to better understand this zoning and better constrain the pressure-temperature-time (PTt) conditions of the metamorphic event of interest.

Major element diffusion in garnet is described by a coupled multicomponent system between four poles: Mg, Fe, Mn and Ca. This can be expressed by a system of parabolic partial differential equations following Fick's first law of diffusion:

```math
\begin{equation}
  \label{ficks_law}
\begin{aligned}
    \frac{\partial C_{i}}{\partial t} = \mathbf{\nabla} \cdot \sum_{j=1}^{n-1} \textbf{D}_{ij} \mathbf{\nabla} C_j,
\end{aligned}
\end{equation}
```

with ``C_i`` the molar fraction of the ``i^{th}`` pole of garnet, ``n`` the total number of poles, and ``\textbf{D}_{ij}`` the interdiffusion coefficient tensor.

For a system of ``n`` components, only ``n-1`` equations are needed, as the sum of the molar fractions of all components is equal to 1. Thus, Ca is made dependent on the other concentrations.

According to Lasaga (1979) [[1]](@ref refs), ``\textbf{D}_{ij}`` can be defined, assuming an ideal system for garnet, with:

```math
\begin{equation}
  \label{diffusion_coefficient}
\begin{aligned}
    \textbf{D}_{ij} = D_i^* \delta_{ij} - \frac{C_i z_i z_j D_i^*}{\sum_{k=1}^{n} C_k z_k^2 D_k^*} (D_j^* - D_{n}^*) =
\begin{bmatrix}
    D_{MgMg} & D_{MgFe} & D_{MgMn} \\
    D_{FeMg} & D_{FeFe} &  D_{FeMn} \\
    D_{MnMg} & D_{MnFe} &  D_{MnMn} \\
\end{bmatrix},
\end{aligned}
\end{equation}
```

with ``D_i^*`` the tracer diffusion coefficient of the ``i^{th}`` component, ``\delta_{ij}`` the Kronecker delta, ``z`` the charge of the species, and ``D_n^*`` the tracer diffusion coefficient of the dependent molar fraction (here Ca).

The tracer diffusion coefficient ``D_i^*`` can be calculated from an Arrhenius relationship:

```math
\begin{equation}
  \label{arrhenius}
D_i^* = D_{0,i} \exp \left( -\frac{E_{a,i} - (P-1)\Delta V^+_i}{RT} \right),
\end{equation}
```

with ``D_{0,i}`` the pre-exponential constant, ``E_{a,i}`` the activation energy of diffusion, ``\Delta V^+_i`` the activation volume of diffusion at 1 bar, ``R`` the universal gas constant, ``T`` the temperature, and ``P`` the pressure.

The pre-exponential constant, ``E_{a,i}``, and ``\Delta V^+_i`` are those of Chakraborty & Ganguly (1992) [[2]](@ref refs). The tracer diffusion coefficient of Ca is defined as ``0.5D_{Fe}``, following the approach of Loomis et al. (1985) [[3]](@ref refs).

### Numerical approach

By defining the PT conditions of the metamorphic event of interest, (3) can be solved for each component, and the diffusion coefficient tensor can be calculated using (2) from the initial major composition data. In DiffusionGarnet.jl, (1) is then discretised in space using finite differences, and the resulting system of ordinary differential equations is solved with ROCK2, a stabilised explicit method (Abdulle & Medovikov (2001) [[4]](@ref refs)) using the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) ecosystem.

## [References](@id refs)

[1] Lasaga, A. C. (1979). Multicomponent exchange and diffusion in silicates. Geochimica et Cosmochimica Acta, 43(4), 455-469.

[2] Chakraborty, S., & Ganguly, J. (1992). Cation diffusion in aluminosilicate garnets: experimental determination in spessartine-almandine diffusion couples, evaluation of effective binary diffusion coefficients, and applications. Contributions to Mineralogy and petrology, 111(1), 74-86.

[3] Loomis, T. P., Ganguly, J., Elphick, S. C., 1985. Experimental determinations of cation diffusitivities in aluminosilicate garnets. II. Multicomponent simulation and tracer diffusion coefficients. Contributions to Mineralogy and Petrology 90, 45–51.

[4] Abdulle, A., & Medovikov, A. A. (2001). Second order Chebyshev methods based on orthogonal polynomials. Numerische Mathematik, 90, 1-18.
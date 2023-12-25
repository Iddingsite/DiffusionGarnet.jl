# [Background](@id background)

Garnet is a mineral commonly used in metamorphic petrology to better understand geological processes, as it occurs in a variety of rock types. This mineral often exhibits a wide range of compositional zoning, which can be interpreted as recording ranges of pressure (P) and temperature (T) conditions. Modelling diffusion processes can help to better understand this zoning and better constrain the pressure-temperature-time (PTt) conditions of the metamorphic event of interest.

Major element diffusion in garnet is described by a coupled multicomponent system between four poles: $Mg$, $Fe$, $Mn$ and $Ca$. This can be expressed by a system of parabolic partial differential equations following Fick's first law of diffusion:

```math
\begin{aligned}
    \frac{\partial C_{i}}{\partial t} = \mathbf{\nabla} \cdot \sum_{j=1}^{n-1} \textbf{D}_{ij} \mathbf{\nabla} C_j
\end{aligned}
```
with $C_i$ the molar fraction of the $i^{th}$ pole of garnet, $n$ the total number of poles, and $\textbf{D}_{ij}$ the interdiffusion coefficient matrix.

For a system of $n$ components, only $n-1$ equations are needed, as the sum of the molar fractions of all components is equal to 1. Thus, $Ca$ is made dependent on the other concentrations.

According to Lasaga (1979) [[1]](@ref refs), $\textbf{D}_{ij}$ can be defined, assuming an ideal system for garnet, with:

```math
\begin{aligned}
    \textbf{D}_{ij} = D_i \delta_{ij} - \frac{C_i D_i}{\sum_{k=1}^{n} C_i D_i} (D_j - D_{n}) =
\begin{bmatrix}
    D_{MgMg} & D_{MgFe} & D_{MgMn} \\
    D_{FeMg} & D_{FeFe} &  D_{FeMn} \\
    D_{MnMg} & D_{MnFe} &  D_{MnMn} \\
\end{bmatrix},
\end{aligned}
```

with $D_i$ the tracer diffusion coefficient of the $i^{th}$ component, $\delta_{ij}$ the Kronecker delta, and $D_n$ the diffusion coefficient of the dependent molar fraction (here $Ca$).

The tracer diffusion coefficient $D_i$ can be calculated with an Arrhenius relationship:

```math
D_i = D_{0,i} \exp \left( -\frac{E_{a,i} - (P-1)\Delta V^+_i}{RT} \right),
```

with $D_{0,i}$ the pre-exponential constant, $E_{a,i}$ the activation energy of diffusion, $\Delta V^+_i$ the activation volume of diffusion at 1 bar, $R$ the universal gas constant, $T$ the temperature, and $P$ the pressure.

The pre-exponential constant $D_{0,i}$, $E_{a,i}$, and $\Delta V^+_i$ are those of Chakraborty and Ganguly (1992) [[2]](@ref refs). The tracer diffusion coefficient of $Ca$ is defined as $0.5D_{Fe}$, following the approach of Loomis et al. (1985) [[3]](@ref refs).

## [References](@id refs)

[1] Lasaga, A. C. (1979). Multicomponent exchange and diffusion in silicates. Geochimica et Cosmochimica Acta, 43(4), 455-469.

[2] Chakraborty, S., & Ganguly, J. (1992). Cation diffusion in aluminosilicate garnets: experimental determination in spessartine-almandine diffusion couples, evaluation of effective binary diffusion coefficients, and applications. Contributions to Mineralogy and petrology, 111(1), 74-86.

[3] Loomis, T. P., Ganguly, J., Elphick, S. C., 1985. Experimental determinations of cation diffusitivities in aluminosilicate garnets. II. Multicomponent simulation and tracer diffusion coefficients. Contributions to Mineralogy and Petrology 90, 45–51.

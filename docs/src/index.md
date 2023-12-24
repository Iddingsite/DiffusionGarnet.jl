# DiffusionGarnet.jl

Garnet is a mineral commonly used in metamorphic petrology to better understand geological processes, as it occurs in a variety of rock types. This mineral often exhibits a wide range of compositional zoning, which can be interpreted as recording ranges of pressure (P) and temperature (T) conditions. Modelling diffusion processes can help to better understand this zoning and better constrain the pressure-temperature-time (PTt) conditions of the metamorphic event of interest.

DiffusionGarnet is a Julia package that models the coupled diffusion of major elements on real garnet data. It currently supports 1D, spherical and 2D coordinates for evenly spaced data and will soon be extended to support 3D coordinates. This package can be used as a teaching tool or for research purposes.

## Installation

DiffusionGarnet may be installed directly from the REPL:
```julia-repl
julia>]
  pkg> add DiffusionGarnet
  pkg> test DiffusionGarnet
```

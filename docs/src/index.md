# DiffusionGarnet.jl

Garnet is a mineral commonly used in metamorphic petrology to better understand geological processes, as it occurs in a variety of different rock types. This mineral often exhibits a wide range of compositional zoning, which could been interpreted as recording ranges of pressure (P) and temperature (T) conditions. Modelling diffusion processes can help to better understand this zoning and better constrain the pressure-temperature-time (PTt) conditions of the metamorphic event of interest.

DiffusionGarnet is a Julia package that can be used to model coupled diffusion of major elements on real garnet data. It currently supports 1D, spherical and 2D coordinates for evenly spaced data and is soon to be extended to support 3D coordinates.

## Installation

DiffusionGarnet may be installed directly from the REPL:
```julia-repl
julia>]
  pkg> add DiffusionGarnet
  pkg> test DiffusionGarnet
```

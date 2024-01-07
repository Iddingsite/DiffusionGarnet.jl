# DiffusionGarnet.jl

This is the documentation page for [DiffusionGarnet.jl](https://github.com/Iddingsite/DiffusionGarnet.jl), a high-level [Julia](https://julialang.org/) package to do coupled diffusion modelling of major elements on real garnet data. It currently supports 1D, spherical, 2D and 3D coordinates for evenly spaced data.

It is built on top of the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package ecosystem and uses [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to allow the user to define appropriate units for their problems. For 2D and 3D models, it uses [ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl) to support multithreading on CPU and parallel computing on GPU.

This package can be used as a teaching tool or for research purposes.

## Installation

DiffusionGarnet.jl is a registered package and may be installed directly from the REPL:

```julia-repl
julia>]
  pkg> add DiffusionGarnet
  pkg> test DiffusionGarnet
```

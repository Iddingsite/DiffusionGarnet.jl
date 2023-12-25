# DiffusionGarnet.jl

This is the documentation page for [DiffusionGarnet.jl](https://github.com/Iddingsite/DiffusionGarnet.jl), a [Julia](https://julialang.org/) package that models coupled diffusion of major elements on real garnet data. It currently supports 1D, spherical and 2D coordinates for evenly spaced data and will soon be extended to support 3D coordinates.

It is built on top of the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package ecosystem and uses [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to allow the user to define appropriate units for their problems. For the 2D models, it uses [ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl) to support multithreading on CPU and parallel computing on GPU.

This package can be used as a teaching tool or for research purposes.

## Installation

DiffusionGarnet.jl is a [registered package](http://pkg.julialang.org) and may be installed directly from the REPL:

```julia-repl
julia>]
  pkg> add DiffusionGarnet
  pkg> test DiffusionGarnet
```

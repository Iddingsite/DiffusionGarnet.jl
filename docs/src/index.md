```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: DiffusionGarnet.jl
  text: Model chemical diffusion in garnet for petrological problems.
  tagline: Made for performance for 1D, 2D, and 3D problems
  actions:
    - theme: brand
      text: Getting started
      link: diffusion_1D
    - theme: alt
      text: Background
      link: background
    - theme: alt
      text: View on Github
      link: https://github.com/Iddingsite/DiffusionGarnet.jl/tree/main
    - theme: alt
      text: API Reference
      link: reference
  image:
    src: logo.png
    alt: DiffusionGarnet.jl

features:
  - icon: 🚀
    title: Fast
    details: Use the package ParallelStencil.jl to allow high performance on both CPU and GPU in 2D and 3D.
    link: /diffusion_3D_GPU

  - icon: 🧩
    title: Composable and Flexible
    details: Built on top of the DifferentialEquations.jl ecosystem, providing flexibility and composability in defining both numerical solvers and problems.

  - icon: 🪨
    title: Compare and Calibrate
    details: Compare easily different diffusion coefficient datasets from the literature with your dataset.
    link: /comparing_diffcoef
---
---
```

# What is DiffusionGarnet.jl?

DiffusionGarnet.jl is a high-level Julia package to model linear or coupled chemical diffusion in trace and major elements using natural garnet data. It is suitable for forward modelling but also inverse modelling to retrieve timescales for geospeedometry. It currently supports 1D, spherical for both evenly and non-evenly spaced data, and 2D and 3D coordinates for evenly spaced data.

## How to install DiffusionGarnet.jl?

DiffusionGarnet.jl is a registered package and may be installed directly with the following command in the Julia REPL

```julia-repl
julia>]
  pkg> add DiffusionGarnet
  pkg> test DiffusionGarnet
```

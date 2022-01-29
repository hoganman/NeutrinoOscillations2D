# NeutrinoOscillations2D

This is a simple Julia script to model two-dimensional flavor oscillations, which is a quantum mechanical phenomenon. The derivation for the ODE solution is given in the included PDF.

## Requirements
 * Julia v1.7

 ```julia
using Pkg
Pkg.add("Plots")
Pkg.add("PyPlots")
```
## Results

The simulation should produce an interactive plot as shown below. 

![Results of numerical simulation](https://github.com/hoganman/NeutrinoOscillations2D/blob/main/FlavorOscillations.png?raw=true)
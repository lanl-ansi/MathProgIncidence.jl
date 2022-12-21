# JuMPIn
**JuMP** **In**cidence. Tools for constructing and analyzing the
incidence graph or matrix of variables and constraints in a JuMP model.

These tools can be used to detect whether and (approximately) why the
Jacobian of equality constraints is structurally or numerically singular,
which commonly happens as the result of a modeling error.

## Installation
TODO. Eventually, I would like the following command to work:
```julia
]add JuMPIn
```

## Dependencies
This package depends on
[JuMP](https://github.com/jump-dev/jump.jl),
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl),
and [BipartiteMatching.jl](https://github.com/IsaacRudich/BipartiteMatching.jl).

## Example
TODO. Create some model with a bug in it, use Dulmage-Mendelsohn to debug.
```julia
using JuMP
import JuMPIn as ji

m = Model()
# TODO: Add constraints and variables to model
con_dmp, var_dmp = ji.dulmage_mendelsohn(m)
oc_con = cat(con_dmp.overconstrained, con_dmp.unmatched, dims = 1)
oc_var = var_dmp.overconstrained

println("Underconstrained subsystem")
println("--------------------------")
println("Constraints")
for con in oc_con
    # TODO: How to print the name of the constraints?
end
println("Variables")
for var in oc_var
    # TODO: How to print the name of the variables?
end
```

## Citation
WIP

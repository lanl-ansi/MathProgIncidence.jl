# MathProgIncidence.jl
**Math Programming Incidence** Graph Analysis.
Tools for constructing and analyzing the
incidence graph or matrix of variables and constraints in a JuMP model.

These tools can be used to detect whether and (approximately) why the
Jacobian of equality constraints is structurally or numerically singular,
which commonly happens as the result of a modeling error.
See the [documentation](https://lanl-ansi.github.io/MathProgIncidence.jl/dev/)
for more information and examples.

## Installation
MathProgIncidence is registered on the Julia general package registry.
It can be installed with:
```julia
julia> ]
(@v1.X) pkg> add MathProgIncidence
```

## Dependencies
This package depends on
[JuMP](https://github.com/jump-dev/jump.jl),
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl),
and the dependencies thereof.

## Example
```julia
using JuMP
import MathProgIncidence as MPIN

m = Model()
comps = [1, 2, 3]
@variable(m, x[comps], start=1/3.0)
@variable(m, flow_comp[comps], start=10.0)
@variable(m, flow, start=30.0)
@variable(m, rho, start=1.0)

@constraint(m, sum_comp_eqn, sum(x) == 1)
@constraint(m, comp_dens_eqn, x*rho .== [1.0, 1.1, 1.2])
@NLconstraint(m, bulk_dens_eqn, 1/rho - sum(1/x[j] for j in comps) == 0)
@constraint(m, comp_flow_eqn, x.*flow .== flow_comp)

igraph = MPIN.IncidenceGraphInterface(m)
con_dmp, var_dmp = MPIN.dulmage_mendelsohn(igraph)
oc_con = [con_dmp.overconstrained..., con_dmp.unmatched...]
oc_var = var_dmp.overconstrained
uc_con = con_dmp.underconstrained
uc_var = [var_dmp.unmatched..., var_dmp.underconstrained...]

println("Overconstrained subsystem")
println("-------------------------")
println("Constraints")
for con in oc_con
    println("  $con")
end
println("Variables")
for var in oc_var
    println("  $var")
end
println()

println("Underconstrained subsystem")
println("--------------------------")
println("Constraints")
for con in uc_con
    println("  $con")
end
println("Variables")
for var in uc_var
    println("  $var")
end
```

## License
MathProgIncidence.jl is open-source software released under the 3-clause BSD license.
See LICENSE.md for more information.

## Citation
If you use MathProgIncidence.jl in your research, we would appreciate you citing
the following paper:
```bibtex
@article{parker2023dulmage,
title = {Applications of the {Dulmage-Mendelsohn} decomposition for debugging nonlinear optimization problems},
journal = {Computers \& Chemical Engineering},
volume = {178},
pages = {108383},
year = {2023},
issn = {0098-1354},
doi = {https://doi.org/10.1016/j.compchemeng.2023.108383},
url = {https://www.sciencedirect.com/science/article/pii/S0098135423002533},
author = {Robert B. Parker and Bethany L. Nicholson and John D. Siirola and Lorenz T. Biegler},
}
```

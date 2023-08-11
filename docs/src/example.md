# Simple Example
This page walks through a simple example of using the Dulmage-Mendelsohn
partition to debug a structural singularity.

## Dulmage-Mendelsohn
We start with some imports and by creating a JuMP model we are interested in.
Usually the model we are interested in debugging is much larger and more
complicated than this. This particular system appeared when debugging a dynamic
1-D partial differential-algebraic equation (PDAE) model representing a chemical
looping combustion reactor.
```julia
using JuMP
import MathProgIncidence

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
```
To check this model for structural singularity, we apply the Dulmage-Mendelsohn
partition.
```julia
igraph = MathProgIncidence.IncidenceGraphInterface(m)
con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
```
If any variables or constraints are unmatched, the (Jacobian of the) model
is structurally singular.
```julia
println("Unmatched constraints")
println("---------------------")
for con in con_dmp.unmatched
    println("  $con")
end
println("Unmatched variables")
println("-------------------")
for var in var_dmp.unmatched
    println("  $var")
end
```
```console
Unmatched constraints
---------------------
(1.0 / rho - (1.0 / x[1] + 1.0 / x[2] + 1.0 / x[3])) - 0.0 = 0
Unmatched variables
-------------------
flow_comp[1]
```
This model has one unmatched constraint and one unmatched variable, so it is
structurally singular. However, the unmatched constraint and variable are not
unique. For example, `flow_comp[2]` could have been unmatched instead of
`flow_comp[1]`.

Unique subsets of variables and constraints that are useful when
debugging a structural singularity are the underconstrained and overconstrained
subsystems. The variables in the underconstrained subsystems are contained in
the `unmatched` and `underconstrained` fields, while the constraints are
contained in the `underconstrained` field.
The variables in the overconstrained subsystem are contained in the
`overconstrained` field, while the constraints are contained in the
`overconstrained` and `unmatched` fields.

We now construct the underconstrained and overconstrained subsystems.
```julia
oc_con = cat(con_dmp.overconstrained, con_dmp.unmatched, dims = 1)
oc_var = var_dmp.overconstrained
uc_con = con_dmp.underconstrained
uc_var = cat(var_dmp.unmatched, var_dmp.underconstrained, dims = 1)
```
And display the constraints and variables contained in each.
```julia
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
```console
Overconstrained subsystem
--------------------------
Constraints
  sum_comp_eqn : x[1] + x[2] + x[3] = 1.0
  comp_dens_eqn : x[1]*rho = 1.0
  comp_dens_eqn : x[2]*rho = 1.1
  comp_dens_eqn : x[3]*rho = 1.2
  (1.0 / rho - (1.0 / x[1] + 1.0 / x[2] + 1.0 / x[3])) - 0.0 = 0
Variables
  x[1]
  rho
  x[2]
  x[3]

Underconstrained subsystem
--------------------------
Constraints
  comp_flow_eqn : x[1]*flow - flow_comp[1] = 0.0
  comp_flow_eqn : x[2]*flow - flow_comp[2] = 0.0
  comp_flow_eqn : x[3]*flow - flow_comp[3] = 0.0
Variables
  flow_comp[1]
  flow
  flow_comp[2]
  flow_comp[3]
```
At this point we must use our intuition about the system being modeled to
identify "what is causing" the singularity.
Looking at the under and over-constrained systems, it appears that we are
missing an equation to calculate `flow`, the total flow rate, and that `rho`,
the bulk density, is over-specified as it is computed by both the bulk density
equation and one of the component density equations.

With this knowledge, we can eventually figure out that (a) we need an equation
to calculate `flow` from density, and that (b) our "bulk density equation" is
actually a *skeletal* density equation. Admittedly, this is difficult to figure
out without the full context behind this particular system.

The following script constructs a new version of the model and checks it for
structural singularity:
```julia
using JuMP
import MathProgIncidence

m = Model()
comps = [1, 2, 3]
@variable(m, x[comps], start=1/3.0)
@variable(m, flow_comp[comps], start=10.0)
@variable(m, flow, start=30.0)
@variable(m, rho_bulk, start=1.0)
@variable(m, rho_skel, start=1.0)
@variable(m, porosity, start=0.25)
velocity = 1.0

@constraint(m, sum_comp_eqn, sum(x) == 1)
@constraint(m, comp_dens_eqn, x*rho_bulk .== [1.0, 1.1, 1.2])
@NLconstraint(m, skel_dens_eqn, 1/rho_skel - sum(1/x[j] for j in comps) == 0)
@constraint(m, bulk_dens_eqn, rho_bulk == (1 - porosity)*rho_skel)
@constraint(m, flow_eqn, flow == velocity * rho_bulk)
@constraint(m, comp_flow_eqn, x.*flow .== flow_comp)

igraph = MathProgIncidence.IncidenceGraphInterface(m)
con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
@assert isempty(con_dmp.unmatched) && isempty(var_dmp.unmatched)
```

# Overview

## What is MathProgIncidence?
MathProgIncidence is a JuMP extension that provides algorithms for analyzing incidence
graphs or matrices defined by JuMP variables and constraints. The standard
graph that is analyzed is a bipartite graph of variables and constraints,
where an edge exists between a variable and constraint if the variable
participates in the constraint.

## Why is MathProgIncidence useful?
In a large modeling/optimization project, especially one involving several
developers, it is fairly easy to make a mistake designing or implementing
an algebraic model. These mistakes commonly cause singularities in the Jacobian
of equality constraints. The algorithms implemented in MathProgIncidence
allow a modeler to identify irreducible subsets of variables and constraints
that are causing singularities, which can be very useful when debugging a
suspected modeling error.

## When should I suspect a modeling error?
One obvious situation is that the solution of optimization problems yields
solutions that do not make sense.
Another is that optimization problems fail to converge in a reasonable amount
of time.
However, optimization problems can fail to converge for a wide variety of
reasons, and the reason for non-convergence is highly dependent on the
algorithm used and optimization problem being solved.
In the case of nonlinear local optimization, symptoms that are often indicative
of modeling errors are large regularization coefficients and large numbers
of restoration iterations in interior point methods.

## What algorithms does MathProgIncidence implement?
1. The **Dulmage-Mendelsohn partition**, which detects subsets of variables and constraints causing a structural singularity. See [`MathProgIncidence.dulmage_mendelsohn`](@ref).
More algorithms may be implemented in the future.

## What models should these be applied to?
These algorithms should be used for applications where singularity of (a
subsystem of) constraints and variables violates an assumption made by the
mathematical/physical model or the solver.
Examples include:
- Index-1 differential-algebraic equations (DAEs)
- Chemical process flowsheets
- Gas pipeline networks
Typically, models of highly nonlinear phenomena distributed over some network
(or time period) are error-prone to implement and could benefit from these
algorithms.
In addition, nonlinear local optimization algorithms such as
[Ipopt](https://github.com/jump-dev/ipopt.jl)
(which are used as subroutines of global and mixed-integer nonlinear
solvers)
assume that the Jacobian of equality constraints is full row rank.
A convenient way to check this is to fix the variables that
you intend to be degrees of freedom, then check the Jacobian of equality
constraints for singularity. If this Jacobian is nonsingular, the full equality
Jacobian is full row rank. Otherwise, this assumption may be violated, and
these algorithms may help deteremine the reason why.

## How can I cite MathProgIncidence?
If you use MathProgIncidence in your research, we would appreciate you citing
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

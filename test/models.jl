#  ___________________________________________________________________________
#
#  JuMPIn.jl: JuMP Incidence Graph Analysis
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

import JuMP

function make_degenerate_flow_model()
    m = JuMP.Model()
    comps = [1, 2, 3]
    @JuMP.variable(m, x[comps], start=1/3.0)
    @JuMP.variable(m, flow_comp[comps], start=10.0)
    @JuMP.variable(m, flow, start=30.0)
    @JuMP.variable(m, rho, start=1.0)

    # sum_component_eqn
    @JuMP.constraint(m, sum_comp_eqn, sum(x) == 1)
    # component_density_eqn
    @JuMP.constraint(m, comp_dens_eqn, x*rho .== [1.0, 1.1, 1.2])
    # density_eqn
    @JuMP.NLconstraint(m, bulk_dens_eqn, 1/rho - sum(1/x[j] for j in comps) == 0)
    # component_flow_eqn
    @JuMP.constraint(m, comp_flow_eqn, x.*flow .== flow_comp)
    return m
end

function make_simple_model()
    m = JuMP.Model()
    JuMP.@variable(m, x[1:3] >= 0, start = 1.0)
    JuMP.@objective(m, Min, x[1]^2 + 2*x[2]^2 + 3*x[3]^2)
    JuMP.@NLconstraint(m, eq1, x[1]^1.5 * x[2]^2 == x[3])
    JuMP.@constraint(m, ineq1, - x[3] <= - 1.0)
    JuMP.@constraint(m, ineq2, x[1] + 2*x[2] >= 4)
    JuMP.@constraint(m, range1, 0 <= x[2] + x[3] <= 10)
    return m
end

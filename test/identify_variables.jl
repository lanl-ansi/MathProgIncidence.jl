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

module TestIdentifyVariables
import JuMP as jmp
import MathOptInterface as moi
using Test: @testset, @test, @test_throws
using JuMPIn: identify_unique_variables

# Local import of JuMP models for testing
include("models.jl")
using .Models: make_degenerate_flow_model


function test_linear()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:sum_comp_eqn])
    # What is happening when we hash a VariableRef? I.e. how does
    # the model get hashed?
    pred_var_set = Set([m[:x][1], m[:x][2], m[:x][3]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_quadratic()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:comp_dens_eqn][1])
    pred_var_set = Set([m[:rho], m[:x][1]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_nonlinear()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:bulk_dens_eqn])
    pred_var_set = Set([m[:rho], m[:x][1], m[:x][2], m[:x][3]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_nonlinear_with_potential_duplicates()
    m = jmp.Model()
    @jmp.variable(m, var[1:2])
    @jmp.NLconstraint(m, con, var[1]^3 + var[1]*var[2] == 1.0)
    variables = identify_unique_variables(con)
    @test(length(variables) == 2)
    @test(Set(variables) == Set([m[:var][1], m[:var][2]]))
end

function test_several_equalities()
    m = make_degenerate_flow_model()
    constraints = [m[:comp_flow_eqn][1], m[:sum_comp_eqn], m[:bulk_dens_eqn]]
    variables = identify_unique_variables(constraints)
    pred_var_set = Set([
        m[:flow], m[:flow_comp][1], m[:x][1], m[:x][2], m[:x][3], m[:rho]
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_several_constraints_with_ineq()
    m = make_degenerate_flow_model()
    @jmp.constraint(m, ineq1, m[:flow_comp][2] >= 1.0)
    @jmp.constraint(m, ineq2, m[:flow_comp][1]^2 + m[:flow_comp][3]^2 <= 1.0)
    constraints = [
        m[:comp_flow_eqn][1],
        m[:sum_comp_eqn],
        m[:bulk_dens_eqn],
        m[:ineq1],
        m[:ineq2],
    ]
    variables = identify_unique_variables(constraints)
    pred_var_set = Set([
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_model()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy)
    @jmp.NLconstraint(m, dummy^3.0 <= 5)
    # Note that include_inequalities=false by default.
    variables = identify_unique_variables(m)
    pred_var_set = Set([
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_model_with_ineq()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy)
    @jmp.NLconstraint(m, dummy^3.0 <= 5)
    variables = identify_unique_variables(m, include_inequality=true)
    pred_var_set = Set([
        m[:dummy],
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_function_with_variable_squared()
    m = jmp.Model()
    @jmp.variable(m, dummy1)
    @jmp.variable(m, dummy2)
    @jmp.constraint(m, dummy_con, dummy1^2 + dummy1*dummy2 == 2.0)
    fcn = moi.get(m, moi.ConstraintFunction(), m[:dummy_con])
    variables = identify_unique_variables(fcn)
    @test(length(variables) == 2)
    @test(Set(variables) == Set([m[:dummy1].index, m[:dummy2].index]))
end

function test_model_bad_constr()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy[1:2])
    @jmp.constraint(m, vectorcon, dummy in moi.Nonnegatives(2))
    @test_throws(
        TypeError,
        identify_unique_variables(m, include_inequality=true),
    )
end

function test_model_bad_constr_no_ineq()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy[1:2])
    @jmp.constraint(m, vectorcon, dummy in moi.Nonnegatives(2))
    # Note that we don't throw an error because we don't attempt to
    # identify variables in "vectorcon" as we don't recognize it as an
    # equality constraint.
    #
    # NOTE: We now do throw an error. I have added a method of
    # set_implies_equality that throws an error if it is provided with a
    # subtype of moi.VectorSet. This is not ~strictly~ correct. However,
    # I have punted for now.
    @test_throws(
        TypeError,
        identify_unique_variables(m, include_inequality=false),
    )
end

function test_fixing_constraint()
    m = jmp.Model()
    @jmp.variable(m, x[1:2])
    jmp.fix(x[1], 1)
    fixing_con = jmp.FixRef(x[1])
    variables = identify_unique_variables(fixing_con)
    @test length(variables) == 1
    @test variables[1] === x[1]
end

function test_inequality_with_bounds()
    m = jmp.Model()
    @jmp.variable(m, x[1:2])
    @jmp.variable(m, 0 <= y[1:2])
    @jmp.constraint(m, x[1] + 2*x[2]^2 == 1)
    variables = identify_unique_variables(m, include_inequality = true)
    pred_var_set = Set([x[1], x[2], y[1], y[2]])
    @test Set(variables) == pred_var_set
end

function test_two_constraints_same_type()
    m = jmp.Model()
    @jmp.variable(m, x[1:3])
    @jmp.constraint(m, eq1, x[1] + x[2] == 2)
    @jmp.constraint(m, eq2, x[2] + 2*x[3] == 3)
    cons = [eq1, eq2]
    vars = identify_unique_variables(cons)
    pred_var_set = Set([x[1], x[2], x[3]])
    @test Set(vars) == pred_var_set
end

function runtests()
    test_linear()
    test_quadratic()
    test_nonlinear()
    test_nonlinear_with_potential_duplicates()
    test_several_equalities()
    test_several_constraints_with_ineq()
    test_model()
    test_model_with_ineq()
    test_model_bad_constr()
    test_model_bad_constr_no_ineq()
    test_function_with_variable_squared()
    test_fixing_constraint()
    test_inequality_with_bounds()
    test_two_constraints_same_type()
    return
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    TestIdentifyVariables.runtests()
end

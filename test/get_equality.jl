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
import MathOptInterface as MOI
using Test: @test, @test_throws, @testset
using JuMPIn: get_equality_constraints, get_inequality_constraints

include("models.jl") # make_degenerate_flow_model, make_simple_model

function get_flow_model_with_inequalities()
    m = make_degenerate_flow_model()
    @JuMP.constraint(m, ineq1, m[:x][1] >= 0)
    @JuMP.constraint(m, ineq2, m[:x][1]^2 + m[:x][2]^2 <= 0.5)
    @JuMP.NLconstraint(m, ineq3, sqrt(m[:x][3]) >= 0.1)
    return m
end

function test_flow_model_with_inequalities()
    m = get_flow_model_with_inequalities()
    eq_cons = get_equality_constraints(m)
    @test length(eq_cons) == 8
    return
end

function test_with_vector_constraint()
    m = JuMP.Model()
    @JuMP.variable(m, var[1:2])
    @JuMP.constraint(m, con, var in MOI.Nonnegatives(2))
    @test_throws(TypeError, eq_cons = get_equality_constraints(m))
    return
end

function test_with_fixed_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:2])
    # These bound constraints are picked up, then discarded as the
    # one-sided interval sets never "imply equality"
    @JuMP.variable(m, 0 <= y <= 1)
    JuMP.fix(x[1], 1.0)
    eq_cons = get_equality_constraints(m)
    @test length(eq_cons) == 1
    return
end

function test_get_inequality_constraints()
    m = make_simple_model()
    ineq_constraints = get_inequality_constraints(m)
    pred_ineq_set = Set([
        JuMP.LowerBoundRef(m[:x][1]),
        JuMP.LowerBoundRef(m[:x][2]),
        JuMP.LowerBoundRef(m[:x][3]),
        m[:ineq1],
        m[:ineq2],
        m[:range1],
    ])
    @test Set(ineq_constraints) == pred_ineq_set
end

@testset "get-equality" begin
    test_flow_model_with_inequalities()
    test_with_vector_constraint()
    test_with_fixed_variables()
    test_get_inequality_constraints()
end

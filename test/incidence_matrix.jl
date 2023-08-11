#  ___________________________________________________________________________
#
#  MathProgIncidence.jl: Math Programming Incidence Graph Analysis
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
import SparseArrays
using Test: @test, @test_throws, @testset

import MathProgIncidence

include("models.jl") # make_degenerate_flow_model

function test_incidence_matrix_from_constraints_and_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1, 2*x[1] + 3*x[2] == 4)
    @JuMP.NLconstraint(m, eq2, 2*x[3]^1.5*x[2] == 1)
    constraints = [eq1, eq2]
    variables = [x[1], x[2], x[3]]
    imat = MathProgIncidence.incidence_matrix(constraints, variables)
    pred_row = [1, 1, 2, 2]
    pred_col = [1, 2, 2, 3]
    pred_val = [1.0, 1.0, 1.0, 1.0]
    m = 2
    n = 3
    pred_mat = SparseArrays.sparse(pred_row, pred_col, pred_val, m, n)
    @test imat == pred_mat
end

function test_incidence_matrix_from_incidence_graph()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1, 2*x[1] + 3*x[2] == 4)
    @JuMP.NLconstraint(m, eq2, 2*x[3]^1.5*x[2] == 1)
    constraints = [eq1, eq2]
    variables = [x[3], x[2], x[1]]
    igraph = MathProgIncidence.IncidenceGraphInterface(constraints, variables)
    imat = MathProgIncidence.incidence_matrix(igraph)
    pred_row = [1, 1, 2, 2]
    pred_col = [3, 2, 2, 1]
    pred_val = [1.0, 1.0, 1.0, 1.0]
    m = 2
    n = 3
    pred_mat = SparseArrays.sparse(pred_row, pred_col, pred_val, m, n)
    @test imat == pred_mat
end

@testset "incidence-matrix" begin
    test_incidence_matrix_from_constraints_and_variables()
    test_incidence_matrix_from_incidence_graph()
end

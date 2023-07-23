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
using JuMPIn: get_bipartite_incidence_graph

include("models.jl") # make_degenerate_flow_model

function test_get_incidence_graph()
    m = make_degenerate_flow_model()
    graph, con_node_map, var_node_map = get_bipartite_incidence_graph(m)
    A, B, E = graph
    edge_set = Set(E)
    @test(length(E) == length(edge_set))
    # To test the incidence graph:
    # - A contains the correct ConstraintRefs in the correct order
    # - B contains the correct VariableRefs in the correct order
    # - For every constraint, and variable in that constraint, the edge
    #   between the variable and constraint appears in E
    # - The total number of edges is equal to the sum of the numbers of
    #   variables in each constraint.
    @test(A == Vector(1:8))
    @test(B == Vector(9:16))
    @test(length(con_node_map) == 8)
    @test(length(var_node_map) == 8)

    con_node_set = Set(values(con_node_map))
    var_node_set = Set(values(var_node_map))
    @test(con_node_set == Set(A))
    @test(var_node_set == Set(B))

    predicted_edges = [
        (m[:sum_comp_eqn], m[:x][1]),
        (m[:sum_comp_eqn], m[:x][2]),
        (m[:sum_comp_eqn], m[:x][3]),
        (m[:comp_dens_eqn][1], m[:rho]),
        (m[:comp_dens_eqn][2], m[:rho]),
        (m[:comp_dens_eqn][3], m[:rho]),
        (m[:comp_dens_eqn][1], m[:x][1]),
        (m[:comp_dens_eqn][2], m[:x][2]),
        (m[:comp_dens_eqn][3], m[:x][3]),
        (m[:comp_flow_eqn][1], m[:flow]),
        (m[:comp_flow_eqn][2], m[:flow]),
        (m[:comp_flow_eqn][3], m[:flow]),
        (m[:comp_flow_eqn][1], m[:x][1]),
        (m[:comp_flow_eqn][2], m[:x][2]),
        (m[:comp_flow_eqn][3], m[:x][3]),
        (m[:comp_flow_eqn][1], m[:flow_comp][1]),
        (m[:comp_flow_eqn][2], m[:flow_comp][2]),
        (m[:comp_flow_eqn][3], m[:flow_comp][3]),
        (m[:bulk_dens_eqn], m[:rho]),
        (m[:bulk_dens_eqn], m[:x][1]),
        (m[:bulk_dens_eqn], m[:x][2]),
        (m[:bulk_dens_eqn], m[:x][3]),
    ]
    @test(length(E) == length(predicted_edges))
    for (con, var) in predicted_edges
        edge = (con_node_map[con], var_node_map[var])
        @test(edge in edge_set)
    end
    return
end

function test_get_incidence_graph_badconstraint()
    m = make_degenerate_flow_model()
    @JuMP.variable(m, var[1:2])
    @JuMP.constraint(m, vectorcon, var in MOI.Nonnegatives(2))
    @test_throws(
        TypeError,
        get_bipartite_incidence_graph(m, include_inequality=true),
    )
    return
end

function test_include_bound_as_inequality()
    m = JuMP.Model()
    @JuMP.variable(m, 0 <= x[1:2])
    @JuMP.constraint(m, eq1, x[1] + 2*x[2] == 1)
    graph, con_node_map, var_node_map = get_bipartite_incidence_graph(
        m, include_inequality = true
    )
    A, B, E = graph
    @test length(A) == 3
    @test length(B) == 2
    @test length(E) == 4
    pred_con_set = Set([eq1, JuMP.LowerBoundRef(x[1]), JuMP.LowerBoundRef(x[2])])
    @test pred_con_set == keys(con_node_map)
    return
end

function test_construct_from_constraints()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1, x[1] + 2*x[2] == 1)
    @JuMP.constraint(m, eq2, x[2]*x[3] == 0.5)
    cons = [eq1, eq2]
    graph, con_node_map, var_node_map = get_bipartite_incidence_graph(cons)
    A, B, E = graph

    @test con_node_map[eq1] == 1
    @test con_node_map[eq2] == 2

    predicted_edges = [(eq1, x[1]), (eq1, x[2]), (eq2, x[2]), (eq2, x[3])]
    @test length(E) == length(predicted_edges)
    edge_set = Set(E)
    for (con, var) in predicted_edges
        @test (con_node_map[con], var_node_map[var]) in edge_set
    end
    return
end

function test_construct_from_constraints_and_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1, x[1] + 2*x[2] == 1)
    @JuMP.constraint(m, eq2, x[2]*x[3] == 0.5)
    cons = [eq2, eq1]
    vars = [x[3], x[2], x[1]]
    graph, con_node_map, var_node_map = get_bipartite_incidence_graph(cons, vars)
    A, B, E = graph

    @test con_node_map[eq1] == 2
    @test con_node_map[eq2] == 1
    @test var_node_map[x[1]] == 5
    @test var_node_map[x[2]] == 4
    @test var_node_map[x[3]] == 3

    predicted_edges = [(eq1, x[1]), (eq1, x[2]), (eq2, x[2]), (eq2, x[3])]
    @test length(E) == length(predicted_edges)
    edge_set = Set(E)
    for (con, var) in predicted_edges
        @test (con_node_map[con], var_node_map[var]) in edge_set
    end
    return
end

@testset "incidence-graph" begin
    test_get_incidence_graph()
    test_get_incidence_graph_badconstraint()
    test_include_bound_as_inequality()
    test_construct_from_constraints()
    test_construct_from_constraints_and_variables()
end

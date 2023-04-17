#  ___________________________________________________________________________
#
#  JuMPIn.jl: JuMP Incidence Graph Analysis
#
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
#  ___________________________________________________________________________

module TestIncidenceGraph

import JuMP as jmp
import MathOptInterface as moi
using Test: @test, @test_throws
using JuMPIn: get_bipartite_incidence_graph

include("models.jl") # Models
using .Models: make_degenerate_flow_model

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
    @jmp.variable(m, var[1:2])
    @jmp.constraint(m, vectorcon, var in moi.Nonnegatives(2))
    @test_throws(
        TypeError,
        get_bipartite_incidence_graph(m, include_inequality=true),
    )
    return
end

function runtests()
    test_get_incidence_graph()
    test_get_incidence_graph_badconstraint()
    return
end

end # module TestIncidenceGraph

if abspath(PROGRAM_FILE) == @__FILE__
    TestIncidenceGraph.runtests()
end

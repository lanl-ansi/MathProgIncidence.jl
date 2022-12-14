module TestIncidenceGraph

import JuMP as jmp
import MathOptInterface as moi

using Test: @test, @test_throws

include("models.jl")
using .models: make_degenerate_flow_model

include("incidence_graph.jl")
using .IncidenceGraph: get_bipartite_incidence_graph

include("get_equality.jl")
using .GetEquality: get_equality_constraints

include("identify_variables.jl")
using .IdentifyVariables: identify_unique_variables

function test_get_incidence_graph()
    m = make_degenerate_flow_model()
    graph, con_node_map, var_node_map = get_bipartite_incidence_graph(m)
    A, B, E = graph
    edge_set = Set(E)
    println(@test(length(E) == length(edge_set)))
    # To test the incidence graph:
    # - A contains the correct ConstraintRefs in the correct order
    # - B contains the correct VariableRefs in the correct order
    # - For every constraint, and variable in that constraint, the edge
    #   between the variable and constraint appears in E
    # - The total number of edges is equal to the sum of the numbers of
    #   variables in each constraint.
    #
    # However, my returned objects are currently useless, as they are indices
    # into arrays that I don't have access to in this scope.
    println(@test(A == Vector(1:8)))
    println(@test(B == Vector(9:16)))
    println(@test(length(con_node_map) == 8))
    println(@test(length(var_node_map) == 8))

    con_node_set = Set(values(con_node_map))
    var_node_set = Set(values(var_node_map))
    println(@test(con_node_set == Set(A)))
    println(@test(var_node_set == Set(B)))

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
end

function test_get_incidence_graph_badconstraint()
    m = make_degenerate_flow_model()
    @jmp.variable(m, var[1:2])
    @jmp.constraint(m, vectorcon, var in moi.Nonnegatives(2))
    @test_throws(
        TypeError,
        get_bipartite_incidence_graph(m, include_inequality=true),
    )
end

function main()
    test_get_incidence_graph()
    println(test_get_incidence_graph_badconstraint())
end

end # module TestIncidenceGraph

if abspath(PROGRAM_FILE) == @__FILE__
    TestInterface.main()
end

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

using Test: @test, @test_throws, @testset
import Graphs
import JuMP
import Ipopt
import MathProgIncidence
import SparseArrays: sparse

include("models.jl") # make_degenerate_flow_model, make_simple_model

function _test_igraph_fields(igraph, constraints, variables)
    @test Set(variables) == keys(igraph._var_node_map)
    @test Set(constraints) == keys(igraph._con_node_map)

    ncon = length(constraints)
    nvar = length(variables)

    var_nodes = Set(values(igraph._var_node_map))
    con_nodes = Set(values(igraph._con_node_map))
    @test con_nodes == Set(1:ncon)
    @test var_nodes == Set(ncon+1:ncon+nvar)

    @test Set(igraph._nodes[1:ncon]) == Set(constraints)
    @test Set(igraph._nodes[ncon+1:ncon+nvar]) == Set(variables)
end


function test_construct_interface(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)

    variables = [
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:flow],
        m[:rho],
    ]
    constraints = [
        m[:sum_comp_eqn],
        m[:comp_dens_eqn][1],
        m[:comp_dens_eqn][2],
        m[:comp_dens_eqn][3],
        m[:bulk_dens_eqn],
        m[:comp_flow_eqn][1],
        m[:comp_flow_eqn][2],
        m[:comp_flow_eqn][3],
    ]
    _test_igraph_fields(igraph, constraints, variables)
    return nothing
end


function test_construct_interface_rectangular(model_function=make_degenerate_flow_model)
    m = model_function()
    @JuMP.constraint(
        m,
        sum_flow_eqn,
        m[:flow] == sum(m[:flow_comp][:]),
    )
    igraph = MathProgIncidence.IncidenceGraphInterface(m)

    variables = [
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:flow],
        m[:rho],
    ]
    constraints = [
        m[:sum_comp_eqn],
        m[:comp_dens_eqn][1],
        m[:comp_dens_eqn][2],
        m[:comp_dens_eqn][3],
        m[:bulk_dens_eqn],
        m[:comp_flow_eqn][1],
        m[:comp_flow_eqn][2],
        m[:comp_flow_eqn][3],
        m[:sum_flow_eqn],
    ]
    _test_igraph_fields(igraph, constraints, variables)
    return nothing
end


function test_get_adjacent_to_linear_constraint(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con = m[:sum_comp_eqn]
    adjacent = MathProgIncidence.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:x][2], m[:x][3]])
    return nothing
end


function test_get_adjacent_to_quadratic_constraint(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con = m[:comp_dens_eqn][1]
    adjacent = MathProgIncidence.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:rho]])
    return nothing
end


function test_get_adjacent_to_nonlinear_constraint(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con = m[:bulk_dens_eqn]
    adjacent = MathProgIncidence.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:x][2], m[:x][3], m[:rho]])
    return nothing
end


function test_get_adjacent_to_variable(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    var = m[:x][2]
    adjacent = MathProgIncidence.get_adjacent(igraph, var)
    incident_cons = [
        m[:sum_comp_eqn],
        m[:bulk_dens_eqn],
        m[:comp_dens_eqn][2],
        m[:comp_flow_eqn][2],
    ]
    @test Set(adjacent) == Set(incident_cons)
    return nothing
end

function test_maximum_matching(model_function=make_degenerate_flow_model)
    m = model_function()
    function _test_matching(matching)
        @test length(matching) == 7
        for (con, var) in matching
            @test typeof(con) <: JuMP.ConstraintRef
            @test typeof(var) <: JuMP.VariableRef
            @test var in Set(MathProgIncidence.get_adjacent(igraph, con))
            @test con in Set(MathProgIncidence.get_adjacent(igraph, var))
        end
        possibly_unmatched_vars = Set([
            m[:flow_comp][1],
            m[:flow_comp][2],
            m[:flow_comp][3],
            m[:flow],
        ])
        possibly_unmatched_cons = Set([
            m[:comp_dens_eqn][1],
            m[:comp_dens_eqn][2],
            m[:comp_dens_eqn][3],
            m[:bulk_dens_eqn],
            m[:sum_comp_eqn],
        ])
        for con in keys(igraph._con_node_map)
            if !(con in keys(matching))
                @test con in possibly_unmatched_cons
            end
        end
        matched_var_set = Set(values(matching))
        for var in keys(igraph._var_node_map)
            if !(var in matched_var_set)
                @test var in possibly_unmatched_vars
            end
        end
    end
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    matching = MathProgIncidence.maximum_matching(igraph)
    _test_matching(matching)
    matching = MathProgIncidence.maximum_matching(m)
    _test_matching(matching)
    return nothing
end

function test_maximum_matching_matrix()
    # Test dense matrix
    matrix = [
        0 1 0;
        0 0 1;
        1 0 0;
    ]
    matching = MathProgIncidence.maximum_matching(matrix)
    @test matching == Dict(1 => 2, 2 => 3, 3 => 1)
    # Test sparse matrix
    matrix = sparse(matrix)
    matching = MathProgIncidence.maximum_matching(matrix)
    @test matching == Dict(1 => 2, 2 => 3, 3 => 1)
    return nothing
end

function test_dulmage_mendelsohn(model_function=make_degenerate_flow_model)
    m = model_function()
    function _test_dm(con_dmp, var_dmp)
        con_undercon = con_dmp.underconstrained
        con_overcon = cat(con_dmp.overconstrained, con_dmp.unmatched, dims=1)
        @test Set(con_undercon) == Set([
            m[:comp_flow_eqn][1], m[:comp_flow_eqn][2], m[:comp_flow_eqn][3]
        ])
        @test Set(con_overcon) == Set([
            m[:comp_dens_eqn][1],
            m[:comp_dens_eqn][2],
            m[:comp_dens_eqn][3],
            m[:bulk_dens_eqn],
            m[:sum_comp_eqn],
        ])
        @test con_dmp.square == []
        var_undercon = cat(var_dmp.underconstrained, var_dmp.unmatched, dims=1)
        var_overcon = var_dmp.overconstrained
        @test Set(var_undercon) == Set([
            m[:flow_comp][1], m[:flow_comp][2], m[:flow_comp][3], m[:flow]
        ])
        @test Set(var_overcon) == Set([m[:x][1], m[:x][2], m[:x][3], m[:rho]])
        @test var_dmp.square == []
    end
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
    _test_dm(con_dmp, var_dmp)
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(m)
    _test_dm(con_dmp, var_dmp)
    return nothing
end

function test_dulmage_mendelsohn_matrix()
    matrix = [
        1 1 0 1;
        0 0 1 0;
        0 0 0 1;
        0 0 0 1;
    ]
    function _test_dm(rowdm, coldm)
        @test rowdm.underconstrained == [1]
        @test rowdm.square == [2]
        @test sort(vcat(rowdm.overconstrained, rowdm.unmatched)) == [3, 4]
        @test coldm.square == [3]
        @test sort(vcat(coldm.underconstrained, coldm.unmatched)) == [1, 2]
        @test coldm.overconstrained == [4]
    end
    rowdm, coldm = MathProgIncidence.dulmage_mendelsohn(matrix)
    _test_dm(rowdm, coldm)
    rowdm, coldm = MathProgIncidence.dulmage_mendelsohn(sparse(matrix))
    _test_dm(rowdm, coldm)
    return nothing
end

function test_overconstrained_due_to_fixed_variable()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:2])
    @JuMP.constraint(m, x[1] + 2*x[2] == 1)
    @JuMP.constraint(m, 3*x[2] - x[2] == 0)
    JuMP.fix(x[1], 3)
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
    @test length(var_dmp.overconstrained) == 2
    @test length(con_dmp.overconstrained) == 2
    @test length(con_dmp.unmatched) == 1
    return
end

function test_overconstrained_due_to_including_bound()
    m = JuMP.Model()
    @JuMP.variable(m, x)
    @JuMP.variable(m, 0.01 <= y)
    @JuMP.constraint(m, 2*x + y == 1)
    @JuMP.NLconstraint(m, x == sqrt(y))
    igraph = MathProgIncidence.IncidenceGraphInterface(m, include_inequality = true)
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
    @test length(var_dmp.overconstrained) == 2
    @test length(con_dmp.overconstrained) == 2
    @test length(con_dmp.unmatched) == 1
    return
end

function test_interface_from_constraints_and_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1,  x[1] + x[2] == 2)
    @JuMP.constraint(m, eq2, x[3]*x[2] == 1.1)
    constraints = [eq1, eq2]
    variables = [x[1], x[3]]
    igraph = MathProgIncidence.IncidenceGraphInterface(constraints, variables)
    _test_igraph_fields(igraph, constraints, variables)
    return
end

function test_matching_from_constraints_and_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1,  x[1] + x[2] == 2)
    @JuMP.constraint(m, eq2, x[3]*x[2] == 1.1)
    constraints = [eq1, eq2]
    variables = [x[1], x[3]]
    matching = MathProgIncidence.maximum_matching(constraints, variables)
    @test length(matching) == 2
    @test matching[eq1] == x[1]
    @test matching[eq2] == x[3]
    return
end

function test_dulmage_mendelsohn_from_constraints_and_variables()
    m = JuMP.Model()
    @JuMP.variable(m, x[1:3])
    @JuMP.constraint(m, eq1,  x[1] + x[2] == 2)
    @JuMP.constraint(m, eq2, x[3]*x[2] == 1.1)
    constraints = [eq1, eq2]
    variables = [x[1], x[3]]
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(constraints, variables)
    @test con_dmp.unmatched == []
    @test con_dmp.underconstrained == []
    @test con_dmp.overconstrained == []
    @test Set(con_dmp.square) == Set(constraints)
    @test var_dmp.unmatched == []
    @test var_dmp.underconstrained == []
    @test var_dmp.overconstrained == []
    @test Set(var_dmp.square) == Set(variables)
    return
end

function test_one_connected_component_igraph(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    cc = MathProgIncidence.connected_components(igraph)
    @test length(cc.var) == 1
    @test length(cc.con) == 1
    @test length(cc.var[1]) == 8
    @test length(cc.con[1]) == 8
    return
end

function test_multiple_connected_components_igraph()
    m = JuMP.Model()
    JuMP.@variable(m, x[1:5] >= 0)
    JuMP.@constraint(m, ineq1, sum(x) >= 3)
    JuMP.@constraint(m, eq1, x[1] + x[3]^2 == 2)
    JuMP.@constraint(m, eq2, x[2] + x[4]^2 == 4)
    JuMP.@constraint(m, eq3, x[5] == 7)
    function _test_cc(concc, varcc)
        predicted_comps = Set(
            [Set([x[1], x[3], eq1]), Set([x[2], x[4], eq2]), Set([x[5], eq3])]
        )
        @test length(varcc) == 3
        @test length(concc) == 3
        for i in 1:3
            comp = Set(cat(varcc[i], concc[i], dims = 1))
            @test comp in predicted_comps
        end
    end
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    cc = MathProgIncidence.connected_components(igraph)
    _test_cc(cc.con, cc.var)
    cc = MathProgIncidence.connected_components(m)
    _test_cc(cc.con, cc.var)
    return
end

function test_one_connected_component_cons_vars(model_function=make_degenerate_flow_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)
    uc_var = [var_dmp.unmatched..., var_dmp.underconstrained...]
    uc_con = con_dmp.underconstrained
    oc_var = var_dmp.overconstrained
    oc_con = [con_dmp.overconstrained..., con_dmp.unmatched...]
    uc_cc = MathProgIncidence.connected_components(uc_con, uc_var)
    oc_cc = MathProgIncidence.connected_components(oc_con, oc_var)
    @test length(uc_cc.con) == 1
    @test length(uc_cc.var) == 1
    @test length(oc_cc.con) == 1
    @test length(oc_cc.var) == 1
    x = m[:x]
    flow_comp = m[:flow_comp]
    flow = m[:flow]
    rho = m[:rho]
    sum_comp_eqn = m[:sum_comp_eqn]
    comp_dens_eqn = m[:comp_dens_eqn]
    bulk_dens_eqn = m[:bulk_dens_eqn]
    comp_flow_eqn = m[:comp_flow_eqn]
    @test Set(uc_cc.con[1]) == Set(comp_flow_eqn)
    @test Set(oc_cc.con[1]) == Set([comp_dens_eqn..., bulk_dens_eqn, sum_comp_eqn])
    @test Set(uc_cc.var[1]) == Set([flow_comp..., flow])
    @test Set(oc_cc.var[1]) == Set([x..., rho])
    return
end

function test_connected_components_matrix()
    matrix = [
        1 0 1;
        1 0 1;
        0 1 0;
    ]
    function _test_cc(rowcc, colcc)
        # Put components in deterministic order
        rowcc = sort(map(sort, rowcc))
        colcc = sort(map(sort, colcc))
        @test rowcc == [[1, 2], [3]]
        @test colcc == [[1, 3], [2]]
    end
    cc = MathProgIncidence.connected_components(matrix)
    _test_cc(cc.con, cc.var)
    cc = MathProgIncidence.connected_components(sparse(matrix))
    _test_cc(cc.con, cc.var)
    return nothing
end

function test_construct_interface_active_inequalities(model_function=make_simple_model)
    m = model_function()
    JuMP.set_optimizer(m, Ipopt.Optimizer)
    JuMP.optimize!(m)

    # Note that this behavior could change if something changes in Ipopt
    # (although this is not likely)

    igraph = MathProgIncidence.IncidenceGraphInterface(
        m, include_active_inequalities = true, tolerance = 1e-6
    )
    constraints = [m[:eq1], m[:ineq1], m[:ineq2]]
    variables = [m[:x][1], m[:x][2], m[:x][3]]
    _test_igraph_fields(igraph, constraints, variables)

    # The default is to use a tolerance of 0.0. With this tolerance, neither
    # of the inequality constraints are active
    igraph = MathProgIncidence.IncidenceGraphInterface(m, include_active_inequalities = true)
    constraints = [m[:eq1]]
    variables = [m[:x][1], m[:x][2], m[:x][3]]
    _test_igraph_fields(igraph, constraints, variables)
    return
end

function test_active_inequalities_no_solution(model_function=make_simple_model)
    m = model_function()
    @test_throws(JuMP.OptimizeNotCalled, igraph = MathProgIncidence.IncidenceGraphInterface(
            m, include_active_inequalities = true, tolerance = 1e-6
        )
    )
    return
end

function test_bad_arguments(model_function=make_simple_model)
    m = model_function()
    @test_throws(ArgumentError, igraph = MathProgIncidence.IncidenceGraphInterface(
            m, include_active_inequalities = true, include_inequality = true
        )
    )
    return
end

function test_block_triangularize(model_function=make_decomposable_model)
    m = model_function()
    igraph = MathProgIncidence.IncidenceGraphInterface(m)
    function _test_blocks(blocks)
        @test length(blocks) == 2
        @test length(blocks[1][1]) == 2
        @test length(blocks[2][1]) == 1
        @test Set(blocks[1][1]) == Set([m[:eq2], m[:eq3]])
        @test Set(blocks[1][2]) == Set([m[:x][1], m[:x][3]])
        @test blocks[2][1] == [m[:eq1]]
        @test blocks[2][2] == [m[:x][2]]
    end
    blocks = MathProgIncidence.block_triangularize(igraph)
    _test_blocks(blocks)
    vars = [m[:x][i] for i in 1:3]
    cons = [m[:eq1], m[:eq2], m[:eq3]]
    blocks = MathProgIncidence.block_triangularize(cons, vars)
    _test_blocks(blocks)
    blocks = MathProgIncidence.block_triangularize(m)
    _test_blocks(blocks)
end

function test_block_triangularize_matrix()
    matrix = [
        1 0 1;
        1 1 1;
        0 1 0;
    ]
    function _test_blocks(block)
        @test blocks[1] == MathProgIncidence.Subsystem(([3], [2]))
        @test Set(blocks[2][1]) == Set([1, 2])
        @test Set(blocks[2][2]) == Set([1, 3])
    end
    blocks = MathProgIncidence.block_triangularize(matrix)
    _test_blocks(blocks)
    blocks = MathProgIncidence.block_triangularize(sparse(matrix))
    _test_blocks(blocks)
    return
end

function test_bfs_from_objective()
    model = make_simple_model_with_ScalarNonlinearFunction()
    objective = JuMP.objective_function(model)
    x = model[:x]
    # Here we exercise the expression interface with no incidence graph and
    # the default depth=1
    tree = MathProgIncidence.limited_bfs(objective)
    nodeset = Set(Any[objective, x...])
    @test nodeset == convert(Set{Any}, Set(tree._nodes))
    predicted_edgeset = Set(Any[
        (objective, x[1]),
        (objective, x[2]),
        (objective, x[3]),
    ])
    sources = map(e -> tree._nodes[e.src], Graphs.edges(tree._dag))
    dests = map(e -> tree._nodes[e.dst], Graphs.edges(tree._dag))
    actual_edgeset = Set{Any}(zip(sources, dests))
    @test actual_edgeset == predicted_edgeset
    return
end

function test_bfs_from_variable()
    # Now we exercise the igraph interface starting from a variable with depth=2
    model = make_simple_model_with_ScalarNonlinearFunction()
    x = model[:x]
    igraph = MathProgIncidence.IncidenceGraphInterface(model; include_inequality = true)
    tree = MathProgIncidence.limited_bfs(igraph, x[1]; depth = 2)
    nodeset = Set(Any[
        x..., model[:eq1], model[:ineq2], JuMP.LowerBoundRef(x[1])
    ])
    @test nodeset == Set{Any}(tree._nodes)
    predicted_edgeset = Set(Any[
        (x[1], model[:eq1]),
        (x[1], model[:ineq2]),
        (x[1], JuMP.LowerBoundRef(x[1])),
        # Note that the particular edges that appear here are not uniquely defined
        # and may change with a different implementation.
        (model[:eq1], x[2]),
        (model[:eq1], x[3]),
    ])
    sources = map(e -> tree._nodes[e.src], Graphs.edges(tree._dag))
    dests = map(e -> tree._nodes[e.dst], Graphs.edges(tree._dag))
    actual_edgeset = Set{Any}(zip(sources, dests))
    @test actual_edgeset == predicted_edgeset
    treestr = sprint(println, tree)
    # on Windows JuMP uses "==" and ">=" instead of "=" and "≥"...
    @test (
        treestr == """
x[1]
├ eq1 : ((x[1] ^ 1.5) * (x[2]²)) - x[3] = 0
│ ├ x[2]
│ └ x[3]
├ ineq2 : x[1] + 2 x[2] ≥ 4
└ x[1] ≥ 0

"""
        || treestr == """
x[1]
├ eq1 : ((x[1] ^ 1.5) * (x[2]²)) - x[3] == 0
│ ├ x[2]
│ └ x[3]
├ ineq2 : x[1] + 2 x[2] >= 4
└ x[1] >= 0

"""
    )
    return
end

@testset "interface" begin
    test_construct_interface()
    test_construct_interface(make_degenerate_flow_model_with_ScalarNonlinearFunction)
    test_construct_interface_rectangular()
    test_construct_interface_rectangular(make_degenerate_flow_model_with_ScalarNonlinearFunction)
    test_get_adjacent_to_linear_constraint()
    test_get_adjacent_to_linear_constraint(make_degenerate_flow_model_with_ScalarNonlinearFunction)
    test_get_adjacent_to_quadratic_constraint()
    test_get_adjacent_to_quadratic_constraint(make_degenerate_flow_model_with_ScalarNonlinearFunction)
    test_get_adjacent_to_nonlinear_constraint()
    test_get_adjacent_to_nonlinear_constraint(make_degenerate_flow_model_with_ScalarNonlinearFunction)
    test_get_adjacent_to_variable()
    test_get_adjacent_to_variable(make_degenerate_flow_model_with_ScalarNonlinearFunction)

    @testset "maximum-matching" begin
        test_maximum_matching()
        test_maximum_matching(make_degenerate_flow_model_with_ScalarNonlinearFunction)
        test_maximum_matching_matrix()
        test_matching_from_constraints_and_variables()
    end

    @testset "dulmage-mendelsohn" begin
        test_dulmage_mendelsohn()
        test_dulmage_mendelsohn(make_degenerate_flow_model_with_ScalarNonlinearFunction)
        test_dulmage_mendelsohn_matrix()
        test_overconstrained_due_to_fixed_variable()
        test_overconstrained_due_to_including_bound()
        test_dulmage_mendelsohn_from_constraints_and_variables()
    end

    @testset "block-triangularize" begin
        test_block_triangularize()
        test_block_triangularize_matrix()
    end

    test_interface_from_constraints_and_variables()

    @testset "connected-components" begin
        test_multiple_connected_components_igraph()
        test_one_connected_component_igraph()
        test_one_connected_component_igraph(make_degenerate_flow_model_with_ScalarNonlinearFunction)
        test_one_connected_component_cons_vars()
        test_one_connected_component_cons_vars(make_degenerate_flow_model_with_ScalarNonlinearFunction)
        test_connected_components_matrix()
    end

    @testset "BFS" begin
        # Testing printed strings can be a bit fragile, so we test the internal
        # IncidenceSubtree data structures. These data are subject to change, so I
        # don't want to write too many tests that are specific to this implementation.
        test_bfs_from_objective()
        test_bfs_from_variable()
    end

    test_construct_interface_active_inequalities()
    test_construct_interface_active_inequalities(make_simple_model_with_ScalarNonlinearFunction)

    test_active_inequalities_no_solution()
    test_active_inequalities_no_solution(make_simple_model_with_ScalarNonlinearFunction)

    test_bad_arguments()
    test_bad_arguments(make_simple_model_with_ScalarNonlinearFunction)
end

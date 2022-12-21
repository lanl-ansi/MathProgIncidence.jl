include("interface.jl") # module Interface
import .Interface as interface

# Note: Do not import from IncidenceGraphInterface here; this will
# lead to the module being defined multiple times, and causes problems
# in the REPL.
#using .Interface: IncidenceGraphInterface

include("models.jl")
using .models: make_degenerate_flow_model
using Test: @test, @test_throws
import JuMP as jmp


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


function test_construct_interface()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)

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


function test_construct_interface_rectangular()
    m = make_degenerate_flow_model()
    @jmp.constraint(
        m,
        sum_flow_eqn,
        m[:flow] == sum(m[:flow_comp][:]),
    )
    igraph = interface.IncidenceGraphInterface(m)

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


function test_get_adjacent_to_linear_constraint()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    con = m[:sum_comp_eqn]
    adjacent = interface.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:x][2], m[:x][3]])
    return nothing
end


function test_get_adjacent_to_quadratic_constraint()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    con = m[:comp_dens_eqn][1]
    adjacent = interface.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:rho]])
    return nothing
end


function test_get_adjacent_to_nonlinear_constraint()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    con = m[:bulk_dens_eqn]
    adjacent = interface.get_adjacent(igraph, con)
    @test Set(adjacent) == Set([m[:x][1], m[:x][2], m[:x][3], m[:rho]])
    return nothing
end


function test_get_adjacent_to_variable()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    var = m[:x][2]
    adjacent = interface.get_adjacent(igraph, var)
    incident_cons = [
        m[:sum_comp_eqn],
        m[:bulk_dens_eqn],
        m[:comp_dens_eqn][2],
        m[:comp_flow_eqn][2],
    ]
    @test Set(adjacent) == Set(incident_cons)
    return nothing
end


function test_maximum_matching()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    matching = interface.maximum_matching(igraph)
    @test length(matching) == 7
    for (con, var) in matching
        @test typeof(con) <: jmp.ConstraintRef
        @test typeof(var) <: jmp.VariableRef
        @test var in Set(interface.get_adjacent(igraph, con))
        @test con in Set(interface.get_adjacent(igraph, var))
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
    return nothing
end


function test_dulmage_mendelsohn()
    m = make_degenerate_flow_model()
    igraph = interface.IncidenceGraphInterface(m)
    con_dmp, var_dmp = interface.dulmage_mendelsohn(igraph)
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
    return nothing
end


function main()
    test_construct_interface()
    test_construct_interface_rectangular()
    test_get_adjacent_to_linear_constraint()
    test_get_adjacent_to_quadratic_constraint()
    test_get_adjacent_to_nonlinear_constraint()
    test_get_adjacent_to_variable()
    test_maximum_matching()
    test_dulmage_mendelsohn()
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

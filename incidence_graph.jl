module IncidenceGraph
import JuMP as jmp
import MathOptInterface as moi

include("get_equality.jl")
using .GetEquality: get_equality_constraints

include("identify_variables.jl")
using .IdentifyVariables: identify_unique_variables

function get_bipartite_incidence_graph(
    model::jmp.Model;
    include_inequality::Bool=false,
)
    if include_inequality
        constraints = jmp.all_constraints(
            model,
            # TODO: Should this be an optional argument to this function?
            include_variable_in_set_constraints=false,
        )   
    else
        constraints = get_equality_constraints(model)
    end
    # Here we get the incidence graph of the constraints, which will by
    # default include all the variables in these constraints.
    # TODO: Should there be an option to include all the variables in the
    # model, even if it results in empty columns?
    return get_bipartite_incidence_graph(constraints)
end

function get_bipartite_incidence_graph(constraints::Vector{jmp.ConstraintRef})
    variables = identify_unique_variables(constraints)
    # We could build up a variable-index map dynamically to get the incidence
    # in a single loop over the constraints, but this is easier to implement.
    return get_bipartite_incidence_graph(constraints, variables)
end

function get_bipartite_incidence_graph(
    constraints::Vector{jmp.ConstraintRef},
    variables::Vector{jmp.VariableRef},
)
    ncon = length(constraints)
    nvar = length(variables)
    # Note the convention we apply: Constraints take the first 1:ncon
    # nodes, while variables take the remainder.
    con_nodes = Vector(1:ncon)
    var_nodes = Vector(ncon+1:ncon+nvar)

    con_node_map = Dict{jmp.ConstraintRef, Int64}(zip(constraints, con_nodes))
    var_node_map = Dict{jmp.VariableRef, Int64}(zip(variables, var_nodes))

    # This could be violated if a variable or constraint appears multiple times.
    # TODO: Fail more gracefully if this happens.
    @assert(ncon == length(con_node_map))
    @assert(nvar == length(var_node_map))

    edges = Vector{Tuple{Int64, Int64}}()
    for con in constraints
        for var in identify_unique_variables(con)
            push!(edges, (con_node_map[con], var_node_map[var]))
        end
    end
    graph = (con_nodes, var_nodes, edges)
    # Return maps as well as the graph so the caller does not need to know
    # about the node convention.
    # While it does not matter for this function, the maps are more useful
    # than vectors of var/con nodes, as the calling functions won't necessarily
    # have ordered vectors of constraints and variables.
    return graph, con_node_map, var_node_map
end

end # module IncidenceGraph

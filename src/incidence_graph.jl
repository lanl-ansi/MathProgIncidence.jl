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

"""
Utility functions for getting the incidence graph of JuMP constraints and
variables.

"""

import JuMP

import JuMPIn: get_equality_constraints, identify_unique_variables

const GraphDataTuple = Tuple{
    # (A, B, E) describing the bipartite graph
    Tuple{Vector{Int}, Vector{Int}, Vector{Tuple{Int, Int}}},
    # Map from constraints to nodes (in A)
    Dict{JuMP.ConstraintRef, Int},
    # Map from variables to nodes (in B)
    Dict{JuMP.VariableRef, Int},
}

"""
    get_bipartite_incidence_graph(model, include_inequality = false)

Return the bipartite incidence graph of (scalar) variables and constraints
in the JuMP model.

The `include_inequality` argument determines whether inequality constraints
(constraints with non-singleton sets) should be included in the graph.

This function returns a tuple `(graph, con_node_map, var_node_map)`
- `graph` -- a tuple `(A, B, E)` where `A` and `B` contain the integer nodes
  in the bipartite sets for constraints and variables and `E` contains
  edges in the form of tuples of integers `(a, b)`, where `a` is in 
  `A` and `b` is in `B`.
- `con_node_map` -- a `Dict` mapping JuMP `ConstraintRef`s to nodes
- `var_node_map` -- a `Dict` mapping JuMP `VariableRef`s to nodes

The constraints in the graph are all the (by default, equality) constraints in
the model, and the variables are those that participate in these constraints.

# Example
```julia-repl
julia> using JuMP

julia> import JuMPIn as ji

julia> m = Model();

julia> @variable(m, v[1:3]);

julia> @constraint(m, eq_1, v[1] + v[3]^2 == 1.0);

julia> @NLconstraint(m, eq_2, v[1]*v[2]^1.5 == 2.0);

julia> graph, con_node_map, var_node_map = ji.get_bipartite_incidence_graph(m);

julia> A, B, E = graph;

julia> M = length(A);

julia> N = length(B);

julia> imat = zeros(M, N);

julia> for (a, b) in E
           imat[a, b-M] = 1.0;
       end

julia> display(imat)
2×3 Matrix{Float64}:
 1.0  1.0  0.0
 0.0  1.0  1.0

```

# Convention
The returned graph follows a convention where, for a JuMP model with
N constraints and M variables, the first N nodes are constraints
and nodes N+1 through N+M are variables. Nodes are always contiguous
integers starting at 1.

# Methods
Methods are implemented that accept a JuMP `Model`, a vector of
`ConstraintRef`s, and vectors of `ConstraintRef`s and `VariableRef`s.
If variables are provided, then only these variables participate in the graph,
regardless of which variables participate in the constraints.

"""
function get_bipartite_incidence_graph(
    model::JuMP.Model;
    include_inequality::Bool = false,
)::GraphDataTuple
    if include_inequality
        # Note that this may generate some constraints that are incompatible
        # with downstream function calls (e.g. constraints involving vector
        # expressions).
        #
        # This is also repeated in identify_unique_variables(Model).
        # Identifying all constraints (including inequalities, but probably
        # not including VectorFunction constraints) may be something we want
        # to standardize at some point. E.g. a get_scalar_constraints function.
        constraints = JuMP.all_constraints(
            model,
            include_variable_in_set_constraints = true,
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

function get_bipartite_incidence_graph(
    constraints::Vector{<:JuMP.ConstraintRef}
)::GraphDataTuple
    variables = identify_unique_variables(constraints)
    # We could build up a variable-index map dynamically to get the incidence
    # in a single loop over the constraints, but this is easier to implement.
    return get_bipartite_incidence_graph(constraints, variables)
end

function get_bipartite_incidence_graph(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::GraphDataTuple
    ncon = length(constraints)
    nvar = length(variables)
    # Note the convention we apply: Constraints take the first 1:ncon
    # nodes, while variables take the remainder.
    con_nodes = Vector(1:ncon)
    var_nodes = Vector(ncon+1:ncon+nvar)

    con_node_map = Dict{JuMP.ConstraintRef, Int64}(zip(constraints, con_nodes))
    var_node_map = Dict{JuMP.VariableRef, Int64}(zip(variables, var_nodes))

    # This could be violated if a variable or constraint appears multiple times.
    # TODO: Fail more gracefully if this happens.
    @assert(ncon == length(con_node_map))
    @assert(nvar == length(var_node_map))

    edges = Vector{Tuple{Int64, Int64}}()
    for con in constraints
        for var in identify_unique_variables(con)
            if var in keys(var_node_map)
                push!(edges, (con_node_map[con], var_node_map[var]))
            end
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

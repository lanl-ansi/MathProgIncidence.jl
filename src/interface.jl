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

"""
A JuMP interface to the algorithms implemented by MathProgIncidence

"""

import JuMP

import MathProgIncidence: get_bipartite_incidence_graph, maximum_matching, GraphDataTuple

import Graphs

"""
Utility function to convert a tuple of nodes and edges into a
Graphs.jl graph
"""
function _tuple_to_graphs_jl(bip_graph)
    A, B, E = bip_graph
    # Assumption here is that A and B are disjoint. And also form
    # a partition of 1:nv.
    nv = length(A) + length(B)
    graph = Graphs.Graph(Graphs.Edge.(E))
    # If E does not cover all vertices, some vertices may not appear in the
    # graph. Add these missing vertices.
    n_missing = nv - Graphs.nv(graph)
    Graphs.add_vertices!(graph, n_missing)
    # Note that by constructing this graph, we have lost our particular
    # bipartition.
    return graph
end

function _maps_to_nodes(con_map, var_map)
    n_nodes = length(con_map) + length(var_map)
    nodes = Vector{Any}([nothing for _ in 1:n_nodes])
    for (con, idx) in con_map
        nodes[idx] = con
    end
    for (var, idx) in var_map
        nodes[idx] = var
    end
    if any(node === nothing for node in nodes)
        throw(Exception)
    end
    return nodes
end

"""
    IncidenceGraphInterface(
        model;
        include_inequality = false,
        include_active_inequalities = false,
        tolerance = 0.0,
    )

A bipartite incidence graph of JuMP constraints and variables.

This is the primary data type accepted by the algorithms implemented in
the remainder of this module.
This type can be instantiated with a JuMP model or a tuple of
`(graph, con_node_map, var_node_map)`, as returned by
[`get_bipartite_incidence_graph`](@ref).
If a model is provided, optional arguments are the same as those provided
to `get_bipartite_incidence_graph`.

Note that the fields of this struct are private, and may change behavior
in a future release without warning.

# Example using only equality constraints
```julia
using JuMP
import MathProgIncidence
m = Model()
@variable(m, v[1:3] >= 0)
@constraint(m, eq_1, v[1] + v[3]^2 == 1.0)
@NLconstraint(m, eq_2, v[1]*v[2]^1.5 == 2.0)
graph = MathProgIncidence.IncidenceGraphInterface(m)
```

# Example including active inequality constraints
```julia
using JuMP
import Ipopt
import MathProgIncidence
m = Model(Ipopt.Optimizer)
@variable(m, v[1:3] >= 0)
@NLconstraint(m, eq_1, v[1]*v[2]^1.5 == 2.0)
@constraint(m, ineq_1, v[1] + v[2] + v[3] >= 7)
@objective(m, Min, v[1]^2 + 2*v[2]^2 + 3*v[3]^2)
optimize!(m)
graph = MathProgIncidence.IncidenceGraphInterface(
    m, include_active_inequalities = true, tolerance = 1e-6
)
```
"""
struct IncidenceGraphInterface
    _graph
    _con_node_map
    _var_node_map
    _nodes
end

IncidenceGraphInterface(
    args::GraphDataTuple
) = IncidenceGraphInterface(
    _tuple_to_graphs_jl(args[1]),
    args[2],
    args[3],
    _maps_to_nodes(args[2], args[3]),
)

IncidenceGraphInterface(
    m::JuMP.Model;
    include_inequality::Bool = false,
    include_active_inequalities::Bool = false,
    tolerance::Float64 = 0.0,
) = IncidenceGraphInterface(
    get_bipartite_incidence_graph(
        m,
        include_inequality = include_inequality,
        include_active_inequalities = include_active_inequalities,
        tolerance = tolerance,
    )
)

IncidenceGraphInterface(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
) = IncidenceGraphInterface(
    get_bipartite_incidence_graph(constraints, variables)
)

"""
    get_adjacent(
        igraph::IncidenceGraphInterface,
        constraint::JuMP.ConstriantRef,
    )::Vector{JuMP.VariableRef}

Return the variables adjacent to a constraint in an incidence graph.

"""
function get_adjacent(
    igraph::IncidenceGraphInterface,
    constraint::JuMP.ConstraintRef,
)::Vector{JuMP.VariableRef}
    con_node = igraph._con_node_map[constraint]
    var_nodes = Graphs.neighbors(igraph._graph, con_node)
    variables = [igraph._nodes[n] for n in var_nodes]
    return variables
end

"""
    get_adjacent(
        igraph::IncidenceGraphInterface,
        variable::JuMP.VariableRef,
    )::Vector{JuMP.ConstraintRef}
    
Return the constraints adjacent to a variable in an incidence graph.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, v[1:3]);

julia> @constraint(m, eq_1, v[1] + v[3] == 1);

julia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);

julia> igraph = MathProgIncidence.IncidenceGraphInterface(m);

julia> adj_cons = MathProgIncidence.get_adjacent(igraph, v[1]);

julia> display(adj_cons)
2-element Vector{ConstraintRef}:
 eq_1 : v[1] + v[3] = 1.0
 v[1] * v[2] ^ 3.0 - 2.0 = 0

```

"""
function get_adjacent(
    igraph::IncidenceGraphInterface,
    variable::JuMP.VariableRef,
)::Vector{JuMP.ConstraintRef}
    var_node = igraph._var_node_map[variable]
    con_nodes = Graphs.neighbors(igraph._graph, var_node)
    constraints = [igraph._nodes[n] for n in con_nodes]
    return constraints
end

"""
    maximum_matching(igraph::IncidenceGraphInterface)::Dict

Compute a maximum matching of variables and constraints in the incidence graph.

The returned `Dict` maps JuMP `ConstraintRef`s to their matched `VariableRef`s.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, v[1:3]);

julia> @constraint(m, eq_1, v[1] + v[3] == 1);

julia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);

julia> igraph = MathProgIncidence.IncidenceGraphInterface(m);

julia> matching = MathProgIncidence.maximum_matching(igraph);

julia> display(matching)
Dict{ConstraintRef, VariableRef} with 2 entries:
  v[1] * v[2] ^ 3.0 - 2.0 = 0 => v[2]
  eq_1 : v[1] + v[3] = 1.0 => v[1]

```

"""
function maximum_matching(
    igraph::IncidenceGraphInterface
)::Dict{JuMP.ConstraintRef, JuMP.VariableRef}
    ncon = length(igraph._con_node_map)
    nodes = igraph._nodes
    con_node_set = Set(1:ncon) # Relying on graph convention here.
    matching = maximum_matching(igraph._graph, con_node_set)
    # matching: constraint nodes -> variable nodes
    jump_matching = Dict(nodes[c] => nodes[v] for (c, v) in matching)
    return jump_matching
end


function maximum_matching(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::Dict{JuMP.ConstraintRef, JuMP.VariableRef}
    igraph = IncidenceGraphInterface(constraints, variables)
    return maximum_matching(igraph)
end

const DMConPartition = NamedTuple{
    (:underconstrained, :square, :overconstrained, :unmatched),
    Tuple{
        Vector{JuMP.ConstraintRef},
        Vector{JuMP.ConstraintRef},
        Vector{JuMP.ConstraintRef},
        Vector{JuMP.ConstraintRef},
    },
}

const DMVarPartition = NamedTuple{
    (:unmatched, :underconstrained, :square, :overconstrained),
    Tuple{
        Vector{JuMP.VariableRef},
        Vector{JuMP.VariableRef},
        Vector{JuMP.VariableRef},
        Vector{JuMP.VariableRef},
    },
}

"""
    dulmage_mendelsohn(igraph::IncidenceGraphInterface)

Return the Dulmage-Mendelsohn partition of variables and constraints
in an incidence graph.

The returned type is a `Tuple` of two `NamedTuple`s, `(con_dmp, var_dmp)`.
These `NamedTuple`s have the following fields:

`con_dmp`:
- `underconstrained` -- The constraints matched with variables that *can
  possibly be* unmatched in a maximum cardinality matching
- `square` -- The constraints that cannot possibly be unmatched in a maximum
  matching
- `overconstrained` -- The constraints that are matched, but *can possibly
  be* unmatched in a maximum matching
- `unmatched` -- The constraints that are not matched in the maximum matching
  that was found

`var_dmp`:
- `unmatched` -- The variables that are not matched in the maximum matching
  that was found
- `underconstrained` -- The variables that *can possibly be* unmatched in a
  maximum matching
- `square` -- The variables that cannot possibly be unmatched in a maximum
  matching
- `overconstrained` -- The variables matched with constraints that *can possibly
  be* unmatched in a maximum cardinality matching

The Dulmage-Mendelsohn partition groups nodes in a bipartite graph into three
unique subsets. In the application to constraints and variables, these may be
thought of as:
- the "overconstrained subsystem", which has more constraints than variables,
- the "underconstrained subsystem", which has more variables than constraints,
- and the "square subsystem", which has the same number of variables as
  constraints

In the `NamedTuple`s returned by this function, the constraints in the
overconstrained subsystem are split into `overconstrained` and `unmatched`,
while the variables in the underconstrained subsystem are split into
`underconstrained` and `unmatched`. This is because it is useful to explicitly
check whether there are any unmatched variables and constraints, and also
useful to recover the maximum matching by `zip`-ing corresponding variables
and constraints.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, v[1:4]);

julia> @constraint(m, eq_1, v[1] + v[3] == 1);

julia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);

julia> @constraint(m, eq_3, v[4]^2 == 3);

julia> igraph = MathProgIncidence.IncidenceGraphInterface(m);

julia> con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph);

julia> # Assert that there are no unmatched constraints

julia> @assert isempty(con_dmp.unmatched);

julia> display(var_dmp.unmatched)
1-element Vector{VariableRef}:
 v[3]

julia> display(var_dmp.underconstrained)
2-element Vector{VariableRef}:
 v[1]
 v[2]

julia> display(con_dmp.underconstrained)
2-element Vector{ConstraintRef}:
 eq_1 : v[1] + v[3] = 1.0
 v[1] * v[2] ^ 3.0 - 2.0 = 0
 
julia> display(var_dmp.square)
1-element Vector{VariableRef}:
 v[4]

julia> display(con_dmp.square)
1-element Vector{ConstraintRef}:
 eq_3 : v[4]² = 3.0

julia> # As there are no unmatched constraints, the overconstrained subsystem is empty
```

"""
function dulmage_mendelsohn(
    igraph::IncidenceGraphInterface
)::Tuple{DMConPartition, DMVarPartition}
    ncon = length(igraph._con_node_map)
    con_node_set = Set(1:ncon)
    con_dmp, var_dmp = dulmage_mendelsohn(igraph._graph, con_node_set)
    # Convert tuples into namedtuples with interpretable names, organized in the
    # order: (underconstrained system, square system, overconstrained system)
    nodes = igraph._nodes
    con_dmp = (
        underconstrained = [nodes[n] for n in con_dmp[3]],
        square = [nodes[n] for n in con_dmp[4]],
        overconstrained = [nodes[n] for n in con_dmp[2]],
        unmatched = [nodes[n] for n in con_dmp[1]],
    )
    var_dmp = (
        unmatched = [nodes[n] for n in var_dmp[1]],
        underconstrained = [nodes[n] for n in var_dmp[2]],
        square = [nodes[n] for n in var_dmp[4]],
        overconstrained = [nodes[n] for n in var_dmp[3]],
    )
    return con_dmp, var_dmp
end

function dulmage_mendelsohn(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::Tuple{DMConPartition, DMVarPartition}
    igraph = IncidenceGraphInterface(constraints, variables)
    return dulmage_mendelsohn(igraph)
end

"""
    connected_components(igraph::IncidenceGraphInterface)

Return the connected components of a bipartite incidence graph of constraints
and variables.

The connected components are returned as two vector-of-vectors, containing
the variables in each connected component and the constraints in each
connected component. Note that the input graph is undirected, so there is no
distinction between strongly and weakly connected components.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, x[1:2] >= 0);

julia> @constraint(m, eq1, x[1] == 1);

julia> @constraint(m, eq2, x[2]^2 == 2);

julia> igraph = MathProgIncidence.IncidenceGraphInterface(m);

julia> con_comps, var_comps = MathProgIncidence.connected_components(igraph);

julia> con_comps
2-element Vector{Vector{ConstraintRef}}:
 [eq1 : x[1] = 1]
 [eq2 : x[2]² = 2]

julia> var_comps
2-element Vector{Vector{VariableRef}}:
 [x[1]]
 [x[2]]

```
"""
function connected_components(
    igraph::IncidenceGraphInterface
)::Tuple{Vector{Vector{JuMP.ConstraintRef}}, Vector{Vector{JuMP.VariableRef}}}
    comps = Graphs.connected_components(igraph._graph)
    ncon = length(igraph._con_node_map)
    nodes = igraph._nodes
    con_node_set = Set(1:ncon)
    con_comps = [[nodes[n] for n in comp if n in con_node_set] for comp in comps]
    var_comps = [[nodes[n] for n in comp if !(n in con_node_set)] for comp in comps]
    return con_comps, var_comps
end

"""
    connected_components(constraints, variables)

Return the connected components of a bipartite incidence graph of constraints
and variables.

The method that accepts constraints and variables directly is convenient for
working with the output of the Dulmage-Mendelsohn partition. It is often used
to decompose and help debug the over and under-constrained subsystems.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, x[1:4] >= 0);

julia> @constraint(m, eq1, x[1] + x[3] == 7);

julia> @constraint(m, eq2, x[2]^2 + x[4]^2 == 1);

julia> igraph = MathProgIncidence.IncidenceGraphInterface(m);

julia> con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph);

julia> uc_con = con_dmp.underconstrained;

julia> uc_var = [var_dmp.unmatched..., var_dmp.underconstrained...];

julia> con_comps, var_comps = MathProgIncidence.connected_components(uc_con, uc_var);

julia> con_comps
2-element Vector{Vector{ConstraintRef}}:
 [eq1 : x[1] + x[3] = 7]
 [eq2 : x[2]² + x[4]² = 1]

julia> var_comps
2-element Vector{Vector{VariableRef}}:
 [x[3], x[1]]
 [x[4], x[2]]

```
"""
function connected_components(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::Tuple{Vector{Vector{JuMP.ConstraintRef}}, Vector{Vector{JuMP.VariableRef}}}
    igraph = IncidenceGraphInterface(constraints, variables)
    return connected_components(igraph)
end

"""
    block_triangularize(igraph::IncidenceGraphInterface)::Vector{Tuple{Vector, Vector}}

Return an ordered partition of constraints and variables that puts the incidence
matrix into block-lower triangular form.

The subsets of the partition correspond to diagonal blocks of the incidence matrix.
These subsets are strongly connected components of a directed version of the bipartite
incidence graph and have the strong Hall property (are irreducible). The
block-triangular ordering of these subsets corresponds to a topological order of the
directed acyclic graph of strongly connected components.

The return type is a vector of tuples of vectors of constraints and variables.

!!! warning
    This is a slightly different return type than the `connected_components` method.
    One or both of these APIs may change in the future for consistency.

# Example
```julia-repl
julia> using JuMP; import MathProgIncidence as MPIN

julia> m = Model(); @variable(m, x[1:3] >= 0);

julia> c1 = @constraint(m, x[1] + x[2]^2 + x[3]^3 == 4);

julia> c2 = @constraint(m, x[1] + x[3] == 1);

julia> c3 = @constraint(m, x[1] * x[3]^0.5 == 3);

julia> igraph = MPIN.IncidenceGraphInterface(m);

julia> blocks = MPIN.block_triangularize(igraph);

julia> length(blocks)
2

julia> for (i, (cb, vb)) in enumerate(blocks)
           println()
           println("Block \$i")
           println("-------")
           for (c, v) in zip(cb, vb) # cb and vb have the same length
               println("\$v, \$c")
           end
       end

Block 1
-------
x[1], (x[1] * (x[3] ^ 0.5)) - 3.0 = 0
x[3], x[1] + x[3] = 1

Block 2
-------
x[2], ((x[2]² + x[1]) + (x[3] ^ 3.0)) - 4.0 = 0

```
"""
function block_triangularize(
    igraph::IncidenceGraphInterface
)::Vector{Tuple{Vector{JuMP.ConstraintRef}, Vector{JuMP.VariableRef}}}
    connodeset = Set(values(igraph._con_node_map))
    ncon = length(igraph._con_node_map)
    nvar = length(igraph._var_node_map)
    matching = maximum_matching(igraph._graph, connodeset)
    nmatch = length(matching)
    if nvar != ncon || nvar != nmatch
        @error (
            "block_triangularize only supports square systems with perfect matchings."
            * "Got nvar=$nvar, ncon=$ncon, n. matched = $nmatch"
        )
    end
    connode_blocks = _block_triangularize(igraph._graph, matching)
    # node_blocks is a vector of tuples of vectors: [([cons], [vars]),...]
    blocks = [
        ([igraph._nodes[n] for n in b], [igraph._nodes[matching[n]] for n in b])
        for b in connode_blocks
    ]
    return blocks
end

"""
    block_triangularize(igraph::IncidenceGraphInterface)::Vector{Tuple{Vector, Vector}}

Return an ordered partition of constraints and variables that puts the incidence
matrix into block-lower triangular form.

This method accepts constraints and variables and is useful for performing the
block-triangular decomposition on the well-constrained subsystem from the
Dulmage-Mendelsohn decomposition.

The return type is a vector of tuples of vectors of constraints and variables.

!!! warning
    This is a slightly different return type than the `connected_components` method.
    One or both of these APIs may change in the future for consistency.

# Example
```julia-repl
julia> using JuMP; import MathProgIncidence as MPIN

julia> m = Model(); @variable(m, x[1:5] >= 0);

julia> c1 = @constraint(m, x[1] + x[2]^2 + x[3]^3 + x[4]^4 +x[5]^5 == 5);

julia> c2 = @constraint(m, x[1] == 1);

julia> c3 = @constraint(m, x[1] * x[3]^0.5 == 3);

julia> c4 = @constraint(m, x[1] - x[2] + x[3] - x[4] == 10);

julia> igraph = MPIN.IncidenceGraphInterface(m);

julia> cdmp, vdmp = MPIN.dulmage_mendelsohn(igraph);

julia> blocks = MPIN.block_triangularize(cdmp.square, vdmp.square);

julia> length(blocks)
2

julia> for (i, (cb, vb)) in enumerate(blocks)
           println("Block \$i")
           println("-------")
           for (c, v) in zip(cb, vb) # cb and vb have the same length
               println("\$v, \$c")
           end
           println()
       end
Block 1
-------
x[1], x[1] = 1

Block 2
-------
x[3], (x[1] * (x[3] ^ 0.5)) - 3.0 = 0

```
"""
function block_triangularize(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::Vector{Tuple{Vector{JuMP.ConstraintRef}, Vector{JuMP.VariableRef}}}
    igraph = IncidenceGraphInterface(constraints, variables)
    return block_triangularize(igraph)
end

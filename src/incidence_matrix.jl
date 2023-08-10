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
import SparseArrays
import Graphs

import MathProgIncidence: get_bipartite_incidence_graph, IncidenceGraphInterface

"""
    incidence_matrix(constraints, variables)::SparseMatrixCSC

Return a `SparseMatrixCSC` with entries corresponding to incidence
between the provided constraints and variables.

Rows correspond to constraints, columns correspond to variables, and
matrix entries correspond to edges in the incidence graph. All matrix
entries have values of `1.0`.

"""
function incidence_matrix(
    constraints::Vector{<:JuMP.ConstraintRef},
    variables::Vector{JuMP.VariableRef},
)::SparseArrays.SparseMatrixCSC
    graph, _, _ = get_bipartite_incidence_graph(constraints, variables)
    A, B, E = graph
    M = length(constraints)
    row = Vector{Int64}()
    col = Vector{Int64}()
    val = Vector{Float64}()
    for (i, j) in E
        # NOTE: Here we exploit the graph's convention. This will need to
        # change if we change how the graph is constructed.
        push!(row, i)
        push!(col, j - M)
        push!(val, 1.0)
    end
    return SparseArrays.sparse(row, col, val)
end

"""
    incidence_matrix(igraph::IncidenceGraphInterface)::SparseMatrixCSC

Return the incidence matrix associated with the provided
`IncidenceGraphInterface`.

"""
function incidence_matrix(
    igraph::IncidenceGraphInterface
)::SparseArrays.SparseMatrixCSC
    row = Vector{Int64}()
    col = Vector{Int64}()
    val = Vector{Float64}()
    M = length(igraph._con_node_map)
    for e in Graphs.edges(igraph._graph)
        i = Graphs.src(e)
        j = Graphs.dst(e)
        if typeof(igraph._nodes[i]) <: JuMP.VariableRef
            i, j = j, i
        end
        # Now we know that j is a variable node
        push!(row, i)
        # NOTE: Here we exploit the graph's convention. This will need to
        # change if we change how the graph is constructed, e.g. we add
        # a node for the objective.
        push!(col, j - M)
        push!(val, 1.0)
    end
    return SparseArrays.sparse(row, col, val)
end

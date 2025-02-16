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

import Graphs

function _get_oriented_projection(graph, matching)
    inv_matching = Dict(b => a for (a, b) in matching)
    nodes = collect(keys(matching))
    @assert Set(nodes) == Set(1:length(matching))
    edgepairs = []
    for e in Graphs.edges(graph)
        u = Graphs.src(e)
        v = Graphs.dst(e)
        @assert u != v # Can't have a self-loop in a bipartite graph
        if u in keys(matching)
            a, b = u, v
        elseif v in keys(matching)
            a, b = v, u
        else
            @error "We should never get here. Contact the MathProgIncidence.jl developers."
        end
        bpair = inv_matching[b]
        if a != bpair
            push!(edgepairs, (a, bpair))
        end
    end
    digraph = Graphs.DiGraph()
    Graphs.add_vertices!(digraph, length(matching))
    for (u, v) in edgepairs
        Graphs.add_edge!(digraph, u, v)
    end
    return digraph
end

function block_triangularize(graph::Graphs.Graph, matching::Dict{Int, Int})
    # TODO: Make sure (a) graph is bipartite and (b) matching is perfect
    # digraph is a new graph. If we ever have to relabel vertices to correspond to
    # Graphs.jl's convention, we'll have to return some mapping to help us get the
    # original vertices back.
    digraph = _get_oriented_projection(graph, matching)
    sccs = Graphs.strongly_connected_components_tarjan(digraph)
    # From the Graphs.jl documentation (https://juliagraphs.org/Graphs.jl/stable/algorithms/connectivity/#Graphs.strongly_connected_components-Tuple{Any}),
    # "returned components will be ordered reverse topologically"
    return sccs
end

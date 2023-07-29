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

import Graphs

import JuMPIn: maximum_matching, _is_valid_bipartition


"""
In this function, matching must contain a key for every matched node.
I.e., Set(keys(matching)) == Set(values(matching))
"""
function _get_projected_digraph(
    graph::Graphs.Graph, nodes::Vector, matching::Dict
)
    # Note that we are constructing a graph in a projected space, and must
    # be aware of whether coordinates are in original or projected spaces.
    orig_proj_map = Dict(n => i for (i, n) in enumerate(nodes))
    n_nodes = Graphs.nv(graph)
    n_nodes_proj = length(nodes)
    matched_nodes = keys(matching)
    node_set = Set(nodes) # Set of nodes in the original space
    edge_set = Set{Tuple{Int64, Int64}}()
    for proj_node in 1:n_nodes_proj
        orig_node = nodes[proj_node]
        if orig_node in matched_nodes
            # In-edges from all (other) neighbors of matched node
            for nbr in Graphs.neighbors(graph, matching[orig_node])
                if nbr != orig_node
                    nbr_proj = orig_proj_map[nbr]
                    push!(edge_set, (nbr_proj, proj_node))
                end
            end
        end
        # TODO: Out edges?
    end
    digraph = Graphs.DiGraph(n_nodes_proj)
    for (n1, n2) in edge_set
        Graphs.add_edge!(digraph, n1, n2)
    end
    return digraph, orig_proj_map
end


function _get_reachable_from(digraph::Graphs.DiGraph, nodes::Vector)
    n_nodes = Graphs.nv(digraph)
    source_set = Set(nodes)
    Graphs.add_vertex!(digraph)
    # Note that root needs to be in this scope so it can be accessed in
    # finally block
    root = n_nodes + 1
    bfs_parents = Vector{Int64}()
    try
        for node in nodes
            Graphs.add_edge!(digraph, root, node)
        end
        bfs_parents = Graphs.bfs_parents(digraph, root)
    finally
        Graphs.rem_vertex!(digraph, root)
    end
    reachable = [
        node for (node, par) in enumerate(bfs_parents)
        if par != 0 && node != root && !(node in source_set)
    ]
    return reachable
end


function dulmage_mendelsohn(graph::Graphs.Graph, set1::Set)
    if !_is_valid_bipartition(graph, set1)
        throw(Exception)
    end
    n_nodes = Graphs.nv(graph)
    set2 = setdiff(Set(1:n_nodes), set1)

    # Compute maximum matching between bipartite sets
    # matching: set1 -> set2
    matching = maximum_matching(graph, set1)
    # inv_matching: set2 -> set1
    inv_matching = Dict(n2 => n1 for (n1, n2) in matching)
    matching = merge(matching, inv_matching)
    matched_set = keys(matching)

    # Compute eight sets of nodes:
    # (a) unmatched nodes in set 1
    # (b) unmatched nodes in set 2
    # (c) nodes in set 1 reachable via alternating path from nodes (a)
    # (d) nodes in set 2 matched with nodes (c)
    # (e) nodes in set 2 reachable via alternating path from nodes (b)
    # (f) nodes in set 1 matched with nodes (e)
    # (g) other nodes in set 1
    # (h) other nodes in set 2

    nodes1 = sort([n for n in set1])
    nodes2 = sort([n for n in set2])

    proj_digraph1, proj_map1 = _get_projected_digraph(graph, nodes1, matching)
    proj_digraph2, proj_map2 = _get_projected_digraph(graph, nodes2, matching)

    unmatched1 = [n for n in nodes1 if !(n in matched_set)]
    unmatched2 = [n for n in nodes2 if !(n in matched_set)]

    # We need the unmatched nodes in the coordinates of the projected graph
    proj_unmatched1 = [proj_map1[n] for n in unmatched1]
    proj_unmatched2 = [proj_map2[n] for n in unmatched2]

    proj_reachable1 = _get_reachable_from(proj_digraph1, proj_unmatched1)
    proj_reachable2 = _get_reachable_from(proj_digraph2, proj_unmatched2)
    reachable1 = [nodes1[n] for n in proj_reachable1]
    reachable2 = [nodes2[n] for n in proj_reachable2]

    # Note that matched_with_reachable1 contains nodes in set 2
    matched_with_reachable1 = [matching[n] for n in reachable1]
    matched_with_reachable2 = [matching[n] for n in reachable2]

    filter = Set(cat(
        unmatched1,
        unmatched2,
        reachable1,
        reachable2,
        matched_with_reachable1,
        matched_with_reachable2;
        dims=1,
    ))
    other1 = [n for n in nodes1 if !(n in filter)]
    other2 = [n for n in nodes2 if !(n in filter)]

    return ((
        unmatched1,
        reachable1,
        matched_with_reachable2,
        other1,
    ), (
        unmatched2,
        reachable2,
        matched_with_reachable1,
        other2,
    ))
end

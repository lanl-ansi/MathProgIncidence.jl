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

const UNMATCHED = nothing
MatchedNodeType{T} = Union{T,typeof(UNMATCHED)}

"""
TODO: This should probably be promoted to Graphs.jl
"""
function _is_valid_bipartition(graph::Graphs.Graph, set1::Set)
    n_nodes = Graphs.nv(graph)
    all_nodes = Set(1:n_nodes)
    if !issubset(set1, all_nodes)
        throw(Exception)
    end
    set2 = setdiff(all_nodes, set1)
    for node in set1
        if !issubset(Graphs.neighbors(graph, node), set2)
            return false
        end
    end
    for node in set2
        if !issubset(Graphs.neighbors(graph, node), set1)
            return false
        end
    end
    return true
end

"""
Determine whether an augmenting path exists and mark distances
so we can compute shortest-length augmenting paths in the DFS.
"""
function _hk_augmenting_bfs!(
    graph::Graphs.AbstractGraph{T},
    set1::Vector{T},
    matching::Dict{T,MatchedNodeType{T}},
    distance::Dict{MatchedNodeType{T},Float64},
)::Bool where {T<:Integer}
    # Initialize queue with the unmatched nodes in set1
    queue = Vector{MatchedNodeType{eltype(graph)}}([
        n for n in set1 if matching[n] == UNMATCHED
    ])

    distance[UNMATCHED] = Inf
    for n in set1
        if matching[n] == UNMATCHED
            distance[n] = 0.0
        else
            distance[n] = Inf
        end
    end

    while !isempty(queue)
        n1 = popfirst!(queue)

        # If n1 is (a) matched or (b) in set1
        if distance[n1] < Inf && n1 != UNMATCHED
            for n2 in Graphs.neighbors(graph, n1)
                # If n2 has not been encountered
                if distance[matching[n2]] == Inf
                    # Give it a distance
                    distance[matching[n2]] = distance[n1] + 1

                    # Note that n2 could be unmatched
                    push!(queue, matching[n2])
                end
            end
        end
    end

    found_augmenting_path = (distance[UNMATCHED] < Inf)
    # The distance to UNMATCHED is the length of the shortest augmenting path
    return found_augmenting_path
end

"""
Compute augmenting paths and update the matching
"""
function _hk_augmenting_dfs!(
    graph::Graphs.AbstractGraph{T},
    root::MatchedNodeType{T},
    matching::Dict{T,MatchedNodeType{T}},
    distance::Dict{MatchedNodeType{T},Float64},
)::Bool where {T<:Integer}
    if root != UNMATCHED
        for n in Graphs.neighbors(graph, root)
            # Traverse edges of the minimum-length alternating path
            if distance[matching[n]] == distance[root] + 1
                if _hk_augmenting_dfs!(graph, matching[n], matching, distance)
                    # If the edge is part of an augmenting path, update the
                    # matching
                    matching[root] = n
                    matching[n] = root
                    return true
                end
            end
        end
        # If we could not find a matched edge that was part of an augmenting
        # path, we need to make sure we don't consider this vertex again
        distance[root] = Inf
        return false
    else
        # Return true to indicate that we are part of an augmenting path
        return true
    end
end

"""
    hopcroft_karp_matching(graph::AbstractGraph)::Dict

Compute a maximum-cardinality matching of a bipartite graph via the
[Hopcroft-Karp algorithm](https://en.wikipedia.org/wiki/Hopcroft-Karp_algorithm).

The return type is a dict mapping nodes to nodes. All matched nodes are included
as keys. For example, if `i` is matched with `j`, `i => j` and `j => i` are both
included in the returned dict.

### Performance

The algorithms runs in O((m + n)n^0.5), where n is the number of vertices and
m is the number of edges. As it does not assume the number of edges is O(n^2),
this algorithm is particularly effective for sparse bipartite graphs.

### Arguments

* `graph`: The bipartite `Graph` for which a maximum matching is computed

### Exceptions

* `ArgumentError`: The provided graph is not bipartite

"""
function hopcroft_karp_matching(graph::Graphs.AbstractGraph{T})::Dict{T,T} where {T<:Integer}
    bmap = Graphs.bipartite_map(graph)
    if length(bmap) != Graphs.nv(graph)
        throw(ArgumentError("Provided graph is not bipartite"))
    end
    set1 = [n for n in Graphs.vertices(graph) if bmap[n] == 1]

    # Initialize "state" that is modified during the algorithm
    matching = Dict{eltype(graph),MatchedNodeType{eltype(graph)}}(
        n => UNMATCHED for n in Graphs.vertices(graph)
    )
    distance = Dict{MatchedNodeType{eltype(graph)},Float64}()

    # BFS to determine whether any augmenting paths exist
    while _hk_augmenting_bfs!(graph, set1, matching, distance)
        for n1 in set1
            if matching[n1] == UNMATCHED
                # DFS to update the matching along a minimum-length
                # augmenting path
                _hk_augmenting_dfs!(graph, n1, matching, distance)
            end
        end
    end
    matching = Dict(i => j for (i, j) in matching if j != UNMATCHED)
    return matching
end

function maximum_matching(graph::Graphs.Graph, set1::Set)
    if !_is_valid_bipartition(graph, set1)
        throw(Exception)
    end
    matching = hopcroft_karp_matching(graph)
    matching = Dict(i => j for (i, j) in matching if i in set1)
    return matching
end

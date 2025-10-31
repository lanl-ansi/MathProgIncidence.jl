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

function _limited_bfs(adjacency_list::Vector{Vector{Int}}, root::Int; depth::Int = 1)
    # This is the queue
    nodes = Int[root]
    idx = 1
    # I use a dict rather than an array here to avoid constructing an O(n)
    # array for every (potentially) small BFS in a (potentially) large graph
    parents = Dict(root => 0)
    # We differ from traditional BFS by imposing a depth limit. We store the index
    # (in the queue) corresponding to the start of the current "level".
    maxdepth = depth
    current_depth = 0
    levelstart = 1
    while idx <= length(nodes)
        node = nodes[idx]
        parent_idx = idx
        idx += 1
        for other in adjacency_list[node]
            if other ∉ keys(parents)
                if parent_idx >= levelstart
                    # When we first encounter the child of a node in the current level,
                    # we advance (a) the depth and (b) the level start index.
                    # The level starts at the child node we are about to add. We do the
                    # checking here so we can break before adding the child if doing
                    # so would violate the depth limit.
                    levelstart = length(nodes) + 1
                    current_depth += 1
                    if current_depth > maxdepth
                        # Break before adding this child to the queue
                        break
                    end
                end
                parents[other] = parent_idx
                push!(nodes, other)
            end
        end
        # Once we exceed the depth limit, we must stop processing children
        # of future nodes.
        if current_depth > maxdepth
            break
        end
    end
    parents = map(n -> parents[n], nodes)
    dag = Graphs.tree(parents) # Returns a DiGraph
    return nodes, dag
end

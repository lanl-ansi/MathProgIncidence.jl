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
import BipartiteMatching as BM


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


function maximum_matching(graph::Graphs.Graph, set1::Set)
    if !_is_valid_bipartition(graph, set1)
        throw(Exception)
    end
    n_nodes = Graphs.nv(graph)
    card1 = length(set1)
    nodes1 = sort([node for node in set1])
    set2 = setdiff(Set(1:n_nodes), set1)
    nodes2 = sort([node for node in set2])
    edge_set = Set((n1, n2) for n1 in nodes1 for n2 in Graphs.neighbors(graph, n1))
    amat = BitArray{2}((r, c) in edge_set for r in nodes1, c in nodes2)
    matching, _ = BM.findmaxcardinalitybipartitematching(amat)
    # Translate row/column coordinates back into nodes of the graph
    graph_matching = Dict(nodes1[r] => nodes2[c] for (r, c) in matching)
    return graph_matching
end

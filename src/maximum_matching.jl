module MaximumMatching

import Graphs as gjl
import BipartiteMatching as bpm


"""
TODO: This should probably be promoted to Graphs.jl
"""
function _is_valid_bipartition(graph::gjl.Graph, set1::Set)
    n_nodes = gjl.nv(graph)
    all_nodes = Set(1:n_nodes)
    if !issubset(set1, all_nodes)
        throw(Exception)
    end
    set2 = setdiff(all_nodes, set1)
    for node in set1
        if !issubset(gjl.neighbors(graph, node), set2)
            return false
        end
    end
    for node in set2
        if !issubset(gjl.neighbors(graph, node), set1)
            return false
        end
    end
    return true
end


function maximum_matching(graph::gjl.Graph, set1::Set)
    if !_is_valid_bipartition(graph, set1)
        throw(Exception)
    end
    n_nodes = gjl.nv(graph)
    card1 = length(set1)
    nodes1 = sort([node for node in set1])
    set2 = setdiff(Set(1:n_nodes), set1)
    nodes2 = sort([node for node in set2])
    edge_set = Set((n1, n2) for n1 in nodes1 for n2 in gjl.neighbors(graph, n1))
    amat = BitArray{2}((r, c) in edge_set for r in nodes1, c in nodes2)
    matching, _ = bpm.findmaxcardinalitybipartitematching(amat)
    # Translate row/column coordinates back into nodes of the graph
    graph_matching = Dict(nodes1[r] => nodes2[c] for (r, c) in matching)
    return graph_matching
end

end

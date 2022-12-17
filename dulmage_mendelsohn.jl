import Graphs as gjl
include("maximum_matching.jl") # _is_valid_bipartition, maximum_matching


"""
In this function, matching must contain a key for every matched node.
I.e., Set(keys(matching)) == Set(values(matching))
"""
function _get_projected_digraph(
    graph::gjl.Graph, nodes::Vector, matching::Dict
)
    # Note that we are constructing a graph in a projected space, and must
    # be aware of whether coordinates are in original or projected spaces.
    orig_proj_map = Dict(n => i for (i, n) in enumerate(nodes))
    n_nodes = gjl.nv(graph)
    n_nodes_proj = length(nodes)
    matched_nodes = keys(matching)
    node_set = Set(nodes) # Set of nodes in the original space
    edge_set = Set{Tuple{Int64}}()
    for proj_node in 1:n_nodes_proj
        orig_node = nodes[proj_node]
        if orig_node in matched_nodes
            # In-edges from all (other) neighbors of matched node
            for nbr in gjl.neighbors(graph, matching[orig_node])
                if nbr != orig_node
                    nbr_proj = orig_proj_map[nbr]
                    push!(edge_set, (nbr_proj, proj_node))
                end
            end
        end
        # TODO: Out edges?
    end
    digraph = gjl.DiGraph(n_proj_nodes)
    for (n1, n2) in edge_set
        gjl.add_edge!(digraph, n1, n2)
    end
    return digraph
end


function _get_reachable_from(digraph::gjl.DiGraph, nodes::Vector)
    n_nodes = gjl.nv(digraph)
    gjl.add_vertex!(digraph)
    root = n_nodes + 1
    for node in nodes
        gjl.add_edge!(digraph, root, node)
    end
    # TODO: get nodes reachable by BFS from root.

    gjl.rem_vertex!(digraph, root)
end


function dulmage_mendelsohn(graph::gjl.Graph, set1::Set)
    if !_is_valid_bipartition(graph, set1)
        throw(Exception)
    end
    n_nodes = gjl.nv(graph)
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

    proj_digraph1 = _get_projected_digraph(graph, nodes1, matching)
    proj_digraph2 = _get_projected_digraph(graph, nodes2, matching)

    unmatched1 = [n for n in nodes1 if !(n in matched_set)]
    unmatched2 = [n for n in nodes2 if !(n in matched_set)]
end

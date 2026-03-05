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
Visualization of the decompositions implemented. These depend on the types
defined in interface.jl
"""

_pad_before_count(label, count, pad) = string(rpad(label, pad), count)

function Base.show(io::IO, igraph::IncidenceGraphInterface)
    # Really, this a graph of _expressions_ and variables...
    rowstr = eltype(keys(igraph._con_node_map)) <: JuMP.ConstraintRef ? "Constraints" : "Rows"
    colstr = eltype(keys(igraph._var_node_map)) <: JuMP.VariableRef ? "Variables" : "Columns"
    nrow = length(igraph._con_node_map)
    ncol = length(igraph._var_node_map)
    pad = max(length(rowstr), length(colstr))
    # TODO: I want to print variables and constraints, but I don't want the display to
    # get unwieldy. We commonly have 10k+ variables and constraints
    println(io, "An Incidence Graph of $rowstr and $colstr")
    println(io, "├ " * _pad_before_count(rowstr * ": ", nrow, pad + 2))
    print(io, "└ " * _pad_before_count(colstr * ": ", ncol, pad + 2))
    return
end

function Base.show(io::IO, dm::DulmageMendelsohnDecomposition)
    uc = Subsystem((dm.con.underconstrained, vcat(dm.var.unmatched, dm.var.underconstrained)))
    oc = Subsystem((vcat(dm.con.overconstrained, dm.con.unmatched), dm.var.overconstrained))
    wc = Subsystem((dm.con.square, dm.var.square))
    print_indented = (subsystem, prefix) -> begin
        lines = split(sprint(print, subsystem), "\n")
        for (i, line) in enumerate(lines)
            if i == 1
                println(io)
            elseif i == length(lines)
                print(io, prefix * line)
            else
                println(io, prefix * line)
            end
        end
        return
    end
    println(io, "Dulmage-Mendelsohn Decomposition")
    print(io, "├ Under-constrained Subsystem")
    print_indented(uc, "│ ")
    println(io)
    print(io, "├ Well-constrained Subsystem")
    print_indented(wc, "│ ")
    println(io)
    print(io, "└ Over-constrained Subsystem")
    print_indented(oc, "  ")
    return
end

function Base.show(io::IO, subsystem::Subsystem)
    rowstr = eltype(subsystem.con) <: JuMP.ConstraintRef ? "Constraints" : "Rows"
    colstr = eltype(subsystem.var) <: JuMP.VariableRef ? "Variables" : "Columns"
    nrow = length(subsystem.con)
    ncol = length(subsystem.var)
    pad = max(length(rowstr), length(colstr))
    println(io, "A Subsystem of $rowstr and $colstr")
    println(io, "├ " * _pad_before_count(rowstr * ": ", nrow, pad + 2))
    print(io, "└ " * _pad_before_count(colstr * ": ", ncol, pad + 2))
    return
end

function Base.show(io::IO, cc::ConnectedComponentDecomposition)
    @assert length(cc.rows) == length(cc.columns)
    n_comps = length(cc.rows)
    print_indented = (subsystem, prefix) -> begin
        lines = split(sprint(print, subsystem), "\n")
        for (i, line) in enumerate(lines)
            if i == 1
                println(io)
            elseif i == length(lines)
                print(io, prefix * line)
            else
                println(io, prefix * line)
            end
        end
        return
    end
    println(io, "Connected Components: $n_comps")
    for (i, (row_comp, col_comp)) in enumerate(zip(cc.rows, cc.columns))
        is_last = i == n_comps
        branch = is_last ? "└" : "├"
        prefix = is_last ? "  " : "│ "
        print(io, "$branch Component $i")
        subsystem = Subsystem((row_comp, col_comp))
        print_indented(subsystem, prefix)
        if !is_last
            println(io)
        end
    end
    return
end

function Base.show(io::IO, subsystems::Vector{Subsystem})
    n_subs = length(subsystems)
    println(io, "A Vector of $n_subs Subsystems")
    for (i, sub) in enumerate(subsystems)
        rowstr = eltype(sub.con) <: JuMP.ConstraintRef ? "Constraints" : "Rows"
        colstr = eltype(sub.var) <: JuMP.VariableRef ? "Variables" : "Columns"
        nrow = length(sub.con)
        ncol = length(sub.var)
        pad = max(length(rowstr), length(colstr))
        is_last = i == n_subs
        branch = is_last ? "└" : "├"
        prefix = is_last ? "  " : "│ "
        println(io, "$branch Subsystem $i")
        println(io, prefix * "├ " * _pad_before_count(rowstr * ": ", nrow, pad + 2))
        if is_last
            print(io, prefix * "└ " * _pad_before_count(colstr * ": ", ncol, pad + 2))
        else
            println(io, prefix * "└ " * _pad_before_count(colstr * ": ", ncol, pad + 2))
        end
    end
    return
end

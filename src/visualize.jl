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

function Base.show(io::IO, dm::DulmageMendelsohnDecomposition)
    nvar = sum(map(length, dm.var))
    ncon = sum(map(length, dm.con))
    #uc_var = vcat(dm.var.underconstrained, dm.var.unmatched)
    #uc_con = dm.con.underconstrained
    #oc_var = dm.var.overconstrained
    #oc_con = vcat(dm.con.overconstrained, dm.con.unmatched)
    uc = Subsystem((dm.con.underconstrained, vcat(dm.var.unmatched, dm.var.underconstrained)))
    oc = Subsystem((vcat(dm.con.overconstrained, dm.con.unmatched), dm.var.overconstrained))
    wc = Subsystem((dm.con.square, dm.var.square))
    msg = "Dulmage-Mendelsohn decomposition with $ncon constraints and $nvar variables"
    nchar = length(msg)
    println(io, msg)
    println(io, repeat("-", nchar))
    println(io, "Under-constrained subsystem: $(length(uc.con)) constraints x $(length(uc.var)) variables")
    println(io, uc)
    println(io, "Over-constrained subsystem:  $(length(oc.con)) constraints x $(length(oc.var)) variables")
    println(io, oc)
    println(io, "Well-constrained subsystem:  $(length(wc.con)) constraints x $(length(wc.var)) variables")
    println(io, wc)
    println(io, repeat("-", nchar))
    return
end

function Base.show(io::IO, subsystem::Subsystem)
    nvar = length(subsystem.var)
    ncon = length(subsystem.con)
    msg = "Subsystem with $nvar variables and $ncon constraints"
    println(io, msg)
    println(io, repeat("-", length(msg)))
    println(io, "Variables:")
    for var in subsystem.var
        println(io, "  $var")
    end
    println(io, "Constraints:")
    for con in subsystem.con
        println(io, "  $con")
    end
    return
end

function Base.show(io::IO, subsystems::Vector{Subsystem})
    n_subs = length(subsystems)
    msg = "Set of $n_subs subsystems"
    println(io, msg)
    println(io, repeat("=", length(msg)))
    for (i, sub) in enumerate(subsystems)
        println()
        nvar = length(sub.var)
        ncon = length(sub.con)
        msg = "Subsystem $i: $nvar variables, $ncon constraints"
        println(io, msg)
        println(io, repeat("-", length(msg)))
        println(io, "Variables:")
        for var in sub.var
            println(io, "  $var")
        end
        println(io, "Constraints:")
        for con in sub.con
            println(io, "  $con")
        end
    end
    println(io, repeat("=", length(msg)))
    return
end

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
    uc_var = vcat(dm.var.underconstrained, dm.var.unmatched)
    uc_con = dm.con.underconstrained
    oc_var = dm.var.overconstrained
    oc_con = vcat(dm.con.overconstrained, dm.con.unmatched)
    msg = "Dulmage-Mendelsohn decomposition with $ncon constraints and $nvar variables"
    nchar = length(msg)
    println(io, msg)
    println(io, repeat("-", nchar))
    println(io, "Under-constrained subsystem: $(length(uc_con)) constraints x $(length(uc_var)) variables")
    println(io, "Over-constrained subsystem:  $(length(oc_con)) constraints x $(length(oc_var)) variables")
    println(io, "Well-constrained subsystem:  $(length(dm.con.square)) constraints x $(length(dm.var.square)) variables")
    println(io, repeat("-", nchar))
    return
end

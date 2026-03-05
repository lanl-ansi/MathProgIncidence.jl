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

using Test
using SparseArrays
import MathProgIncidence as MPIN

include("models.jl")

@testset "visualize" begin
    m = make_degenerate_flow_model()
    igraph = MPIN.IncidenceGraphInterface(m)
    igraph_expected = """An Incidence Graph of Constraints and Variables
├ Constraints: 8
└ Variables:   8
"""
    @test sprint(println, igraph) == igraph_expected

    dm = MPIN.dulmage_mendelsohn(igraph)
    dm_expected = """Dulmage-Mendelsohn Decomposition
├ Under-constrained Subsystem
│ ├ Constraints: 3
│ └ Variables:   4
├ Well-constrained Subsystem
│ ├ Rows:    0
│ └ Columns: 0
└ Over-constrained Subsystem
  ├ Constraints: 5
  └ Variables:   4
"""
    @test sprint(println, dm) == dm_expected

    m_decomp = make_decomposable_model()
    blocks = MPIN.block_triangularize(m_decomp)
    blocks_expected = """A Vector of 2 Subsystems
├ Subsystem 1
│ ├ Constraints: 2
│ └ Variables:   2
└ Subsystem 2
  ├ Constraints: 1
  └ Variables:   1
"""
    @test sprint(println, blocks) == blocks_expected

    matrix = sparse([
        1 0 1;
        1 0 1;
        0 1 0;
    ])
    cc = MPIN.connected_components(matrix)
    cc_expected = """Connected Components: 2
├ Component 1
│ ├ Rows:    2
│ └ Columns: 2
└ Component 2
  ├ Rows:    1
  └ Columns: 1
"""
    @test sprint(println, cc) == cc_expected
end

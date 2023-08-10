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

using Documenter
using MathProgIncidence

makedocs(
    sitename = "JuMPIn",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "Overview" => "overview.md",
        "Simple Example" => "example.md",
        "API Reference" => [
            "reference/get_equality.md",
            "reference/identify_variables.md",
            "reference/incidence_graph.md",
            "reference/interface.md",
            "reference/incidence_matrix.md",
        ],
    ],
)

deploydocs(repo = "github.com/lanl-ansi/JuMPIn.jl.git")

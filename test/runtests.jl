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

using JuMPIn
using Test

@testset "JuMPIn.jl" begin

    @testset "IdentifyVariables" begin
        include("identify_variables.jl")
        TestIdentifyVariables.runtests()
    end

    @testset "GetEquality" begin
        include("get_equality.jl")
        TestGetEquality.runtests()
    end

    @testset "IncidenceGraph" begin
        include("incidence_graph.jl")
        TestIncidenceGraph.runtests()
    end

    @testset "Interface" begin
        include("interface.jl")
        TestInterface.runtests()
    end

end

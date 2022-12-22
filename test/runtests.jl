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

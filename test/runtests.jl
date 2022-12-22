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

end

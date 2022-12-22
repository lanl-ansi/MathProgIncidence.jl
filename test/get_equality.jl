module TestGetEquality
import JuMP as jmp
import MathOptInterface as moi
using Test: @test, @test_throws
using JuMPIn: get_equality_constraints

include("models.jl") # Models

function get_flow_model_with_inequalities()
    m = Models.make_degenerate_flow_model()
    # TODO: Some assertions with the @assert macro
    @jmp.constraint(m, ineq1, m[:x][1] >= 0)
    @jmp.constraint(m, ineq2, m[:x][1]^2 + m[:x][2]^2 <= 0.5)
    @jmp.NLconstraint(m, ineq3, sqrt(m[:x][3]) >= 0.1)
    return m
end

function test_flow_model_with_inequalities()
    m = get_flow_model_with_inequalities()
    eq_cons = get_equality_constraints(m)
    println("N. eq con: ", length(eq_cons))
    @test length(eq_cons) == 8
    return
end

function test_with_vector_constraint()
    m = jmp.Model()
    @jmp.variable(m, var[1:2])
    @jmp.constraint(m, con, var in moi.Nonnegatives(2))
    @test_throws(TypeError, eq_cons = get_equality_constraints(m))
end

function runtests()
    test_flow_model_with_inequalities()
    test_with_vector_constraint()
end

end # module TestGetEquality

if abspath(PROGRAM_FILE) == @__FILE__
    TestGetEquality.runtests()
end

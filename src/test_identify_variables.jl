module TestIdentifyVariables
import JuMP as jmp
import MathOptInterface as moi
using Test: @testset, @test, @test_throws
include("models.jl")
using .models: make_degenerate_flow_model
include("identify_variables.jl")
using .IdentifyVariables: identify_unique_variables


function test_linear()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:sum_comp_eqn])
    # What is happening when we hash a VariableRef? I.e. how does
    # the model get hashed?
    pred_var_set = Set([m[:x][1], m[:x][2], m[:x][3]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_quadratic()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:comp_dens_eqn][1])
    pred_var_set = Set([m[:rho], m[:x][1]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_nonlinear()
    m = make_degenerate_flow_model()
    variables = identify_unique_variables(m[:bulk_dens_eqn])
    pred_var_set = Set([m[:rho], m[:x][1], m[:x][2], m[:x][3]])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_nonlinear_with_potential_duplicates()
    m = jmp.Model()
    @jmp.variable(m, var[1:2])
    @jmp.NLconstraint(m, con, var[1]^3 + var[1]*var[2] == 1.0)
    variables = identify_unique_variables(con)
    @test(length(variables) == 2)
    @test(Set(variables) == Set([m[:var][1], m[:var][2]]))
end

function test_several_equalities()
    m = make_degenerate_flow_model()
    constraints = [m[:comp_flow_eqn][1], m[:sum_comp_eqn], m[:bulk_dens_eqn]]
    variables = identify_unique_variables(constraints)
    pred_var_set = Set([
        m[:flow], m[:flow_comp][1], m[:x][1], m[:x][2], m[:x][3], m[:rho]
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_several_constraints_with_ineq()
    m = make_degenerate_flow_model()
    @jmp.constraint(m, ineq1, m[:flow_comp][2] >= 1.0)
    @jmp.constraint(m, ineq2, m[:flow_comp][1]^2 + m[:flow_comp][3]^2 <= 1.0)
    constraints = [
        m[:comp_flow_eqn][1],
        m[:sum_comp_eqn],
        m[:bulk_dens_eqn],
        m[:ineq1],
        m[:ineq2],
    ]
    variables = identify_unique_variables(constraints)
    pred_var_set = Set([
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_model()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy)
    @jmp.NLconstraint(m, dummy^3.0 <= 5)
    # Note that include_inequalities=false by default.
    variables = identify_unique_variables(m)
    pred_var_set = Set([
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_model_with_ineq()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy)
    @jmp.NLconstraint(m, dummy^3.0 <= 5)
    variables = identify_unique_variables(m, include_inequality=true)
    pred_var_set = Set([
        m[:dummy],
        m[:flow],
        m[:flow_comp][1],
        m[:flow_comp][2],
        m[:flow_comp][3],
        m[:x][1],
        m[:x][2],
        m[:x][3],
        m[:rho],
    ])
    @test(length(variables) == length(pred_var_set))
    @test(pred_var_set == Set(variables))
end

function test_function_with_variable_squared()
    m = jmp.Model()
    @jmp.variable(m, dummy1)
    @jmp.variable(m, dummy2)
    @jmp.constraint(m, dummy_con, dummy1^2 + dummy1*dummy2 == 2.0)
    fcn = moi.get(m, moi.ConstraintFunction(), m[:dummy_con])
    variables = identify_unique_variables(fcn)
    @test(length(variables) == 2)
    @test(Set(variables) == Set([m[:dummy1].index, m[:dummy2].index]))
end

function test_model_bad_constr()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy[1:2])
    @jmp.constraint(m, vectorcon, dummy in moi.Nonnegatives(2))
    @test_throws(
        TypeError,
        identify_unique_variables(m, include_inequality=true),
    )
end

function test_model_bad_constr_no_ineq()
    m = make_degenerate_flow_model()
    @jmp.variable(m, dummy[1:2])
    @jmp.constraint(m, vectorcon, dummy in moi.Nonnegatives(2))
    # Note that we don't throw an error because we don't attempt to
    # identify variables in "vectorcon" as we don't recognize it as an
    # equality constraint.
    #
    # NOTE: We now do throw an error. I have added a method of
    # set_implies_equality that throws an error if it is provided with a
    # subtype of moi.VectorSet. This is not ~strictly~ correct. However,
    # I have punted for now.
    @test_throws(
        TypeError,
        identify_unique_variables(m, include_inequality=false),
    )
end

function main()
    println(test_linear())
    println(test_quadratic())
    println(test_nonlinear())
    println(test_nonlinear_with_potential_duplicates())
    println(test_several_equalities())
    println(test_several_constraints_with_ineq())
    println(test_model())
    println(test_model_with_ineq())
    println(test_model_bad_constr())
    println(test_model_bad_constr_no_ineq())
    println(test_function_with_variable_squared())
    return
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    println(test_model_bad_constr_no_ineq())
    #TestIdentifyVariables.main()
end

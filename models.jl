module models
export make_degenerate_flow_model
import JuMP as jmp

function make_degenerate_flow_model()
    m = jmp.Model()
    comps = [1, 2, 3]
    @jmp.variable(m, x[comps], start=1/3.0)
    @jmp.variable(m, flow_comp[comps], start=10.0)
    @jmp.variable(m, flow, start=30.0)
    @jmp.variable(m, rho, start=1.0)

    # sum_component_eqn
    @jmp.constraint(m, sum_comp_eqn, sum(x) == 1)
    # component_density_eqn
    @jmp.constraint(m, comp_dens_eqn, x*rho .== [1.0, 1.1, 1.2])
    # density_eqn
    @jmp.NLconstraint(m, bulk_dens_eqn, 1/rho - sum(1/x[j] for j in comps) == 0)
    # component_flow_eqn
    @jmp.constraint(m, comp_flow_eqn, x.*flow .== flow_comp)
    return m
end

function main()
    # The purpose of this function is so that, to "test" this module's
    # functionality from the command line, I only have to remember
    # to call the function main, rather than any of the module-specific
    # functions.
    println(models.make_degenerate_flow_model())
end

end # module models

if abspath(PROGRAM_FILE) == @__FILE__
    models.main()
end

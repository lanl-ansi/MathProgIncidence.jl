module GetEquality
import JuMP as jmp
import MathOptInterface as moi

function get_set_of_constraint(
    model::jmp.Model,
    constraint::jmp.ConstraintRef,
    index::moi.ConstraintIndex,
)
    set = moi.get(model, moi.ConstraintSet(), constraint)
    return set
end

function get_set_of_constraint(
    model::jmp.Model,
    constraint::jmp.ConstraintRef,
    index::moi.Nonlinear.ConstraintIndex,
)
    nl_con = model.nlp_model.constraints[index]
    return nl_con.set
end

function set_implies_equality(set::moi.EqualTo)::Bool
    return true
end

# Q: How does Julia determine if an argument has the right type.
# I.e. why is this interpreted as "set must be a subtype of Union(...)"
# rather than "set must be an 'instance' of Union()".
#function set_implies_equality(
#    set::Union{moi.GreaterThan, moi.LessThan},
#    # Note that ^this works, but this doesn't
#    # set::Type{<:Union{moi.GreaterThan, moi.LessThan}},
#)::Bool
#    return false
#end

function set_implies_equality(set::T)::Bool where T<:moi.AbstractVectorSet
    throw(TypeError(set, moi.AbstractScalarSet, typeof(set)))
end

# NOTE: Have not implemented tolerance arg in any calling function yet.
function set_implies_equality(
    set::moi.Interval,
    tolerance::Float64=Float64(0)
)::Bool
    return abs(set.upper - set.lower) <= tolerance
end

# NOTE: Overload this function to define other constraint sets (e.g. Zeros)
# as equalities.
#
# I intend this to catch anything that is a subtype of ConstraintSet.
# Is the "subtype syntax" necessary?
# It appears, e.g., moi.GreaterThan is not a subtype of moi.ConstraintSet.
# This is non-intuitive for me...
# So this won't catch arbitrary moi Sets...
# The right type to use is actually AbstractSet?
# Nope. This doesn't catch GreatherThan/LessThan
# ^ Wrong. Actually it does, but my syntax was wrong. I guess providing
# Type{<:SomeType} enforces that the argument must be a Type, rather than
# an instance of the type. "issubclass", where what I'm doing below is
# "isinstance"
function set_implies_equality(
    set::T
)::Bool where T<:moi.AbstractSet
    return false
end

function is_equality(constraint::jmp.ConstraintRef)::Bool
    model = constraint.model
    index = constraint.index
    set = get_set_of_constraint(model, constraint, index)
    return set_implies_equality(set)
end

function get_equality_constraints(
    constraints::Vector{jmp.ConstraintRef}
)::Vector{jmp.ConstraintRef}
    eq_cons = Vector{jmp.ConstraintRef}()
    for con in constraints
        if is_equality(con)
            push!(eq_cons, con)
        end
    end
    return eq_cons
end

function get_equality_constraints(model::jmp.Model)::Vector{jmp.ConstraintRef}
    constraints = jmp.all_constraints(
        model,
        include_variable_in_set_constraints=false,
    )
    return get_equality_constraints(constraints)
end

end # module get_equality

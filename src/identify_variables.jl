module IdentifyVariables
import JuMP as jmp
import MathOptInterface as moi

include("get_equality.jl")
using .GetEquality: get_equality_constraints


# TODO: This file implements functions that filter duplicates from the
# vectors of identified variables. It may be useful at somepoint to
# identify variables, preserving duplicates. This can be implemented
# when/if there is a need.


"""
    identify_unique_variables(
        model::Model, include_inequality::Bool=false
    )

Returns a vector of variables that participate in constraints in the model.

Each variable appears at most one time in the returned vector.
"""
function identify_unique_variables(
    model::jmp.Model; include_inequality::Bool=false,
)::Vector{jmp.VariableRef}
    if include_inequality
        # Note that this may include some constraints which are not compatible
        # with downstream function calls (e.g. constraints where the function
        # and terms are vectors). We will allow downstream errors to be raised
        # rather than silently ignoring these constraints.
        constraints = jmp.all_constraints(
            model,
            # TODO: Should this be an optional argument to this function?
            include_variable_in_set_constraints=false,
        )   
    else
        constraints = get_equality_constraints(model)
    end
    return identify_unique_variables(constraints)
end


"""
    identify_unique_variables(
        constraints::Vector,
    )::Vector{VariableRef}

Returns a vector of variables that participate in the provided constraints.

Each variable appears at most one time in the returned vector.
If we receive a vector, we just assume it is a vector of constraints.
This is because I couldn't get an argument of type Vector{ConstraintRef}
to work...
"""
function identify_unique_variables(
    constraints::Vector,
    # FIXME: Couldn't get this working with Vector{ConstraintRef}...
)::Vector{jmp.VariableRef}
    variables = Vector{jmp.VariableRef}()
    for con in constraints
        append!(variables, identify_unique_variables(con))
    end
    return _filter_duplicates(variables)
end


"""
    identify_unique_variables(constraint::ConstraintRef)

Returns a vector containing the variables that participate in this constraint.

Each variable appears at most one time in the returned vector.
"""
function identify_unique_variables(
    constraint::jmp.ConstraintRef,
)::Vector{jmp.VariableRef}
    return identify_unique_variables(constraint, constraint.index)
end


"""
    identify_unique_variables(
        constraint::ConstraintRef,
        index::ConstraintIndex,
    )::Vector{VariableRef}

Returns a vector containing the variables that participate in this constraint.

Each variable appears at most one time in the returned vector.
The `index` argument is provided so we know whether we have a "regular"
constraint or a nonlinear constraint.
"""
function identify_unique_variables(
    constraint::jmp.ConstraintRef,
    index::moi.ConstraintIndex,
)::Vector{jmp.VariableRef}
    model = constraint.model
    fcn = moi.get(model, moi.ConstraintFunction(), constraint)
    varidxs = identify_unique_variables(fcn)
    varrefs = [jmp.VariableRef(model, idx) for idx in varidxs]
    return varrefs
end


"""
    identify_unique_variables(
        constraint::ConstraintRef,
        index::MOI.Nonlinear.ConstraintIndex,
    )::Vector{VariableRef}

Returns a vector containing the variables that participate in this constraint.

Each variable appears at most one time in the returned vector.
The `index` argument is provided so we know whether we have a "regular"
constraint or a nonlinear constraint.
"""
function identify_unique_variables(
    constraint::jmp.ConstraintRef,
    index::moi.Nonlinear.ConstraintIndex,
)::Vector{jmp.VariableRef}
    model = constraint.model
    nlmod = model.nlp_model
    nlcons = nlmod.constraints
    con = nlcons[index]
    expr = con.expression

    # Identify variables in expression
    # This could be its own function, but I inlined it here to allow a more
    # useful error message if we ever encounter NODE_VARIABLE.
    # Will be moved when/if the need arises.
    variables = Vector{moi.VariableIndex}()
    for node in expr.nodes
        if node.type == moi.Nonlinear.NODE_MOI_VARIABLE
            push!(variables, moi.VariableIndex(node.index))
        elseif node.type == moi.Nonlinear.NODE_VARIABLE
            # I do not know under what situation this could occur, but
            # am throwing this error to be defensive.
            throw(DomainError(
                node.type,
                """Encountered a NODE_VARIABLE while parsing constraint\
                   $constraint,\nbut we do not have an NLPEvaluator\
                   to resolve this into a NODE_MOI_VARIABLE.
                   Something has gone wrong.""",
            ))
        end
    end

    # Return VariableRefs
    refs = [jmp.VariableRef(model, idx) for idx in variables]
    return _filter_duplicates(refs)
end


"""
    identify_unique_variables(
        fcn::Union{ScalarQuadraticFunction, ScalarAffineFunction},
    )::Vector{VariableIndex}

Returns the variables that appear in the provided function.

NOTE: Only ScalarQuadraticFunction and ScalarAffineFunction are supported.
This can be changed there is demand for other functions. For each type of
supported function, the _get_variable_terms function should be defined.
Then, for the type of each term, an additional identify_unique_variables
function should be implemented.
"""
function identify_unique_variables(
    fcn::Union{moi.ScalarQuadraticFunction, moi.ScalarAffineFunction},
)::Vector{moi.VariableIndex}
    variables = Vector{moi.VariableIndex}()
    for terms in _get_variable_terms(fcn)
        for term in terms
            for var in identify_unique_variables(term)
                push!(variables, var)
            end
        end
    end
    return _filter_duplicates(variables)
end

function identify_unique_variables(
    fcn::T
)::Vector{moi.VariableIndex} where {T<:moi.AbstractVectorFunction}
    throw(TypeError(
        fcn,
        Union{moi.ScalarQuadraticFunction, moi.ScalarAffineFunction},
        typeof(fcn),
    ))
end


"""
    _get_variable_terms(fcn)

Returns a tuple of vectors of terms for the provided MOI function.

Currently implemented only for ScalarQuadraticFunction and ScalarAffineFunction.
"""
function _get_variable_terms(fcn::moi.ScalarQuadraticFunction)
    return (fcn.quadratic_terms, fcn.affine_terms)
end

function _get_variable_terms(fcn::moi.ScalarAffineFunction)
    return (fcn.terms,)
end


"""
    identify_unique_variables(term)

Returns the variables that participate in the provided term of a MOI function.

A variable appears at most one time in the returned vector. Currently
implemented only for ScalarAffineTerm and ScalarQuadraticTerm.
"""
function identify_unique_variables(
    term::moi.ScalarAffineTerm,
)::Vector{moi.VariableIndex}
    return [term.variable]
end

function identify_unique_variables(
    term::moi.ScalarQuadraticTerm,
)::Vector{moi.VariableIndex}
    # Note that this compares indices only. If somehow these indices refer to
    # different models, this could be incorrect.
    if term.variable_1 === term.variable_2
        return [term.variable_1]
    else
        return [term.variable_1, term.variable_2]
    end
end


"""
    _filter_duplicates

Returns a vector of variables of indices that does not contain duplicates.
"""
function _filter_duplicates(
    indices::Vector{moi.VariableIndex},
)::Vector{moi.VariableIndex}
    # Note that by hashing only the variable index, we implicitly assume that
    # all variables come from the same model. I believe this is safe.
    seen = Set{moi.VariableIndex}()
    filtered = Vector{moi.VariableIndex}()
    for idx in indices
        if !(idx in seen)
            push!(seen, idx)
            push!(filtered, idx)
        end
    end
    return filtered
end


function _filter_duplicates(
    variables::Vector{jmp.VariableRef},
)::Vector{jmp.VariableRef}
    # What happens when VariableRef gets hashed? I.e. how is the
    # hash of the model computed?
    seen = Set{jmp.VariableRef}()
    filtered = Vector{jmp.VariableRef}()
    for var in variables
        if !(var in seen)
            push!(seen, var)
            push!(filtered, var)
        end
    end
    return filtered
end

end

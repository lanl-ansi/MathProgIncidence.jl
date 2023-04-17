#  ___________________________________________________________________________
#
#  JuMPIn.jl: JuMP Incidence Graph Analysis
#
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
#  ___________________________________________________________________________

"""
Utility functions for identifying variables that participate in constraints.

"""
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
    identify_unique_variables(constraints::Vector)::Vector{JuMP.VariableRef}

Return a vector of variables that participate in the provided constraints.

Each variable appears at most one time in the returned vector.
If we receive a vector, we just assume it is a vector of constraints.
This is because I couldn't get an argument of type `Vector{ConstraintRef}`
to work...

Note that this function is also accessible via the JuMPIn module.

# Example
```julia-repl
julia> using JuMP
julia> import JuMPIn as ji
julia> m = Model();
julia> @variable(m, v[1:3]);
julia> @constraint(m, eq_1, v[2] == 1);
julia> @NLconstraint(m, eq_2, v[2]*v[3]^1.5 == 2);
julia> vars = ji.identify_unique_variables([eq_1, eq_2]);
julia> display(vars)
2-element Vector{VariableRef}:
 v[2]
 v[3]
```

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
    identify_unique_variables(model::JuMP.Model, include_inequality::Bool=false)

Return a vector of variables that participate in constraints in the model.

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
    identify_unique_variables(constraint::JuMP.ConstraintRef)

Return a vector containing the variables that participate in this constraint.

Each variable appears at most one time in the returned vector.

"""
function identify_unique_variables(
    constraint::jmp.ConstraintRef,
)::Vector{jmp.VariableRef}
    return identify_unique_variables(constraint, constraint.index)
end


"""
    identify_unique_variables(
        constraint::JuMP.ConstraintRef,
        index::JuMP.ConstraintIndex,
    )::Vector{JuMP.VariableRef}

Return a vector containing the variables that participate in this constraint.

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
        constraint::JuMP.ConstraintRef,
        index::MathOptInterface.Nonlinear.ConstraintIndex,
    )::Vector{JuMP.VariableRef}

Return a vector containing the variables that participate in this constraint.

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
    identify_unique_variables(fcn)::Vector{JuMP.VariableIndex}

Return the variables that appear in the provided MathOptInterface function.

# Implementation
Only `ScalarQuadraticFunction` and `ScalarAffineFunction` are supported.
This can be changed there is demand for other functions. For each type of
supported function, the `_get_variable_terms` function should be defined.
Then, for the type of each term, an additional `identify_unique_variables`
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

Return a tuple of vectors of terms for the provided MathOptInterface function.

# Implementation
Currently implemented only for `ScalarQuadraticFunction` and
`ScalarAffineFunction`.

"""
function _get_variable_terms(fcn::moi.ScalarQuadraticFunction)
    return (fcn.quadratic_terms, fcn.affine_terms)
end

function _get_variable_terms(fcn::moi.ScalarAffineFunction)
    return (fcn.terms,)
end


"""
    identify_unique_variables(term)

Return the variables that participate in the provided term of a
MathOptInterface function.

A variable appears at most one time in the returned vector.

# Implementation
Currently implemented only for `ScalarAffineTerm` and `ScalarQuadraticTerm`.

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

Return a vector of variables of indices that does not contain duplicates.

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

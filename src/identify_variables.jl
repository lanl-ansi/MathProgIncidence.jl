#  ___________________________________________________________________________
#
#  MathProgIncidence.jl: Math Programming Incidence Graph Analysis
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

"""
Utility functions for identifying variables that participate in constraints.

"""

import JuMP
import MathOptInterface as MOI

import MathProgIncidence: get_equality_constraints


"""
    identify_unique_variables(constraints::Vector{JuMP.ConstraintRef})::Vector{JuMP.VariableRef}

Return a vector of variables that participate in the provided constraints.

Each variable appears at most one time in the returned vector.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, v[1:3]);

julia> @constraint(m, eq_1, v[2] == 1);

julia> @NLconstraint(m, eq_2, v[2]*v[3]^1.5 == 2);

julia> vars = MathProgIncidence.identify_unique_variables([eq_1, eq_2]);

julia> display(vars)
2-element Vector{VariableRef}:
 v[2]
 v[3]

```

"""
function identify_unique_variables(
    constraints::Vector{<:JuMP.ConstraintRef},
)::Vector{JuMP.VariableRef}
    variables = Vector{JuMP.VariableRef}()
    for con in constraints
        append!(variables, identify_unique_variables(con))
    end
    return _filter_duplicates(variables)
end


"""
    identify_unique_variables(
        model::JuMP.Model,
        include_inequality = false,
        include_active_inequalities = false,
        tolerance = 0.0,
    )

Return a vector of variables that participate in constraints in the model.

Each variable appears at most one time in the returned vector.

# Optional keyword arguments
- `include_inequality`: Whether to include variables that participate in *any*
  inequality constraints
- `include_active_inequalities`: Whether to include variables that participate
  *active* inequality constraints (at the latest solution).
  `include_active_inequalities` and `include_inequality` are mutually exclusive.
- `tolerance`: Tolerance used to determine if an inequality constraint is active

"""
function identify_unique_variables(
    model::JuMP.Model;
    include_inequality::Bool=false,
    include_active_inequalities::Bool = false,
    tolerance::Float64 = 0.0,
)::Vector{JuMP.VariableRef}
    # Note that this method exists mostly for convenience, and is not used
    # by any of the "upstream" methods, e.g. get_bipartite_incidence_graph
    constraints = _get_constraints(
        model,
        include_inequality = include_inequality,
        include_active_inequalities = include_active_inequalities,
        tolerance = tolerance,
    )
    return identify_unique_variables(constraints)
end


"""
    identify_unique_variables(constraint::JuMP.ConstraintRef)

Return a vector containing the variables that participate in this constraint.

Each variable appears at most one time in the returned vector.

"""
function identify_unique_variables(
    constraint::JuMP.ConstraintRef,
)::Vector{JuMP.VariableRef}
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
    constraint::JuMP.ConstraintRef,
    index::MOI.ConstraintIndex,
)::Vector{JuMP.VariableRef}
    model = constraint.model
    fcn = MOI.get(model, MOI.ConstraintFunction(), constraint)
    varidxs = identify_unique_variables(fcn)
    varrefs = [JuMP.VariableRef(model, idx) for idx in varidxs]
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
    constraint::JuMP.ConstraintRef,
    index::MOI.Nonlinear.ConstraintIndex,
)::Vector{JuMP.VariableRef}
    model = constraint.model
    nlmod = model.nlp_model
    nlcons = nlmod.constraints
    con = nlcons[index]
    expr = con.expression

    # Identify variables in expression
    # This could be its own function, but I inlined it here to allow a more
    # useful error message if we ever encounter NODE_VARIABLE.
    # Will be moved when/if the need arises.
    variables = Vector{MOI.VariableIndex}()
    for node in expr.nodes
        if node.type == MOI.Nonlinear.NODE_MOI_VARIABLE
            push!(variables, MOI.VariableIndex(node.index))
        elseif node.type == MOI.Nonlinear.NODE_VARIABLE
            # I do not know under what situation this could occur, but
            # am throwing this error to be defensive.
            throw(DomainError(
                node.type,
                """Encountered a NODE_VARIABLE while parsing constraint
                   $constraint,\nbut we do not have an NLPEvaluator
                   to resolve this into a NODE_MOI_VARIABLE.
                   Something has gone wrong.""",
            ))
        end
    end

    # Return VariableRefs
    refs = [JuMP.VariableRef(model, idx) for idx in variables]
    return _filter_duplicates(refs)
end

"""
    identify_unique_variables(fcn)::Vector{JuMP.VariableIndex}

Return the variables that appear in the provided MathOptInterface function.

No variable will appear more than once.

# Implementation
Only `ScalarNonlinearFunction`, `ScalarQuadraticFunction`, and
`ScalarAffineFunction` are supported. This can be changed there is demand for
other functions.
These methods are implemented by first identifying all variables that
participate in the function, then filtering out duplicate variables.

"""
function identify_unique_variables(
    fcn::Union{
        MOI.ScalarAffineFunction,
        MOI.ScalarQuadraticFunction,
        MOI.ScalarNonlinearFunction,
    },
)::Vector{MOI.VariableIndex}
    variables = _identify_variables(fcn)
    return _filter_duplicates(variables)
end

# This method is used to handle function-in-set constraints.
function identify_unique_variables(
    var::MOI.VariableIndex
)::Vector{MOI.VariableIndex}
    return [var]
end

# NOTE: We will get a MethodError if this is called with a non-vector,
# non-Affine/Quadratic/Nonlinear function. Should probably implement
# some default method to catch that.
function identify_unique_variables(
    fcn::MOI.AbstractVectorFunction
)::Vector{MOI.VariableIndex}
    throw(TypeError(
        fcn,
        Union{
            MOI.ScalarAffineFunction,
            MOI.ScalarQuadraticFunction,
            MOI.ScalarNonlinearFunction,
        },
        typeof(fcn),
    ))
end

"""
    _identify_variables(fcn)::Vector{JuMP.VariableIndex}

Return all variables that appear in the provided MathOptInterface function.

Duplicates may be present in the resulting vector of variables.
If there is a use case, these methods can be made public.

# Implementation
Implemented for `ScalarAffineFunction`, `ScalarQuadraticFunction`, and
`ScalarNonlinearFunction`. Relies on an underlying _collect_variables!
method that (potentially recursively) builds up a vector of variables
in-place.

"""
function _identify_variables(
    fcn::Union{
        MOI.ScalarAffineFunction,
        MOI.ScalarQuadraticFunction,
        MOI.ScalarNonlinearFunction,
    },
)::Vector{MOI.VariableIndex}
    variables = Vector{MOI.VariableIndex}()
    _collect_variables!(variables, fcn)
    return variables
end

"""
    _collect_variables!(variables, fcn)::Vector{JuMP.VariableIndex}

Add variables from `fcn` to the `variables` vector

# Implementation
Implemented for `ScalarAffineFunction`, `ScalarQuadraticFunction`, and
`ScalarNonlinearFunction`. For affine and quadratic functions, we iterate
over terms with the `_get_variable_terms` function. For nonlinear functions,
we recurse until pushing a root node onto the variable stack (or hitting
an affine/quadratic function). Methods may need to be added if more node types
are added to `ScalarNonlinearFunction` in the future.

"""
function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    fcn::MOI.ScalarNonlinearFunction,
)::Vector{MOI.VariableIndex}
    for arg in fcn.args
        _collect_variables!(variables, arg)
    end
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    fcn::Union{MOI.ScalarQuadraticFunction, MOI.ScalarAffineFunction},
)::Vector{MOI.VariableIndex}
    for terms in _get_variable_terms(fcn)
        for term in terms
            _collect_variables!(variables, term)
        end
    end
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    var::MOI.VariableIndex,
)::Vector{MOI.VariableIndex}
    push!(variables, var)
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    var::Float64,
)::Vector{MOI.VariableIndex}
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    var::Int,
)::Vector{MOI.VariableIndex}
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    term::MOI.ScalarAffineTerm,
)::Vector{MOI.VariableIndex}
    push!(variables, term.variable)
    return variables
end

function _collect_variables!(
    variables::Vector{MOI.VariableIndex},
    term::MOI.ScalarQuadraticTerm,
)::Vector{MOI.VariableIndex}
    push!(variables, term.variable_1, term.variable_2)
    return variables
end

"""
    _get_variable_terms(fcn)

Return a tuple of vectors of terms for the provided MathOptInterface function.

# Implementation
Currently implemented only for `ScalarQuadraticFunction` and
`ScalarAffineFunction`.

"""
function _get_variable_terms(fcn::MOI.ScalarQuadraticFunction)
    return (fcn.quadratic_terms, fcn.affine_terms)
end

function _get_variable_terms(fcn::MOI.ScalarAffineFunction)
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
    term::MOI.ScalarAffineTerm,
)::Vector{MOI.VariableIndex}
    return [term.variable]
end

function identify_unique_variables(
    term::MOI.ScalarQuadraticTerm,
)::Vector{MOI.VariableIndex}
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
    indices::Vector{MOI.VariableIndex},
)::Vector{MOI.VariableIndex}
    # Note that by hashing only the variable index, we implicitly assume that
    # all variables come from the same model. I believe this is safe.
    seen = Set{MOI.VariableIndex}()
    filtered = Vector{MOI.VariableIndex}()
    for idx in indices
        if !(idx in seen)
            push!(seen, idx)
            push!(filtered, idx)
        end
    end
    return filtered
end

function _filter_duplicates(
    variables::Vector{JuMP.VariableRef},
)::Vector{JuMP.VariableRef}
    # What happens when VariableRef gets hashed? I.e. how is the
    # hash of the model computed?
    seen = Set{JuMP.VariableRef}()
    filtered = Vector{JuMP.VariableRef}()
    for var in variables
        if !(var in seen)
            push!(seen, var)
            push!(filtered, var)
        end
    end
    return filtered
end

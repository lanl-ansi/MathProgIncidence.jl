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
Utility functions for identifying JuMP constraints that define equalities.

"""

import JuMP
import MathOptInterface as MOI

function _get_set_of_constraint(
    model::JuMP.Model,
    constraint::JuMP.ConstraintRef,
    index::MOI.ConstraintIndex,
)
    set = MOI.get(model, MOI.ConstraintSet(), constraint)
    return set
end

function _get_set_of_constraint(
    model::JuMP.Model,
    constraint::JuMP.ConstraintRef,
    index::MOI.Nonlinear.ConstraintIndex,
)
    nl_con = model.nlp_model.constraints[index]
    return nl_con.set
end

function set_implies_equality(set::MOI.EqualTo)::Bool
    return true
end

function set_implies_equality(set::MOI.AbstractVectorSet)::Bool
    throw(TypeError(set, MOI.AbstractScalarSet, typeof(set)))
end

# NOTE: Have not implemented tolerance arg in any calling function yet.
function set_implies_equality(
    set::MOI.Interval,
    tolerance::Float64=Float64(0)
)::Bool
    return abs(set.upper - set.lower) <= tolerance
end

"""
    set_implies_equality(set::MOI.AbstractSet)

Detect whether the set defines an equality constraint, i.e. is a singleton.

# Implementation
Methods are defined for the following `MathOptInterface.Set`s:
- `MathOptInterface.EqualTo` 
- `MathOptInterface.Interval`

If a `MathOptInterface.AbstractVectorSet` is provided, an error is raised.
For any other type of set, `false` is returned. To support additional
types of constraints in [`is_equality`](@ref) and
[`get_equality_constraints`](@ref), additional methods of
`set_implies_equality` should be defined.

"""
function set_implies_equality(set::MOI.AbstractSet)::Bool
    return false
end

function set_implies_inequality(set::MOI.GreaterThan)::Bool
    return true
end

function set_implies_inequality(set::MOI.LessThan)::Bool
    return true
end

function set_implies_inequality(
    set::MOI.Interval;
    tolerance::Float64=0.0,
)::Bool
    # NOTE: tolerance has not been implemented in any calling function
    return abs(set.upper - set.lower) > tolerance
end

"""
    set_implies_inequality(set::MOI.AbstractSet)

Detect whether the set defines an inequality constraint.

This function is defined for scalar sets. Calling with a vector set
will result in a `TypeError`.
"""
function set_implies_inequality(set::MOI.AbstractSet)::Bool
    return false
end

function set_implies_inequality(set::MOI.AbstractVectorSet)::Bool
    throw(TypeError(set, MOI.AbstractScalarSet, typeof(set)))
end

# TODO: set_implies_inequality(set::MOI.AbstractSet, tolerance)
# that just calls set_implies_inequality(set). I.e. the default
# is to ignore the tolerance.

"""
    is_equality(constraint::JuMP.ConstraintRef)::Bool

Detect whether a constraint is an equality constraint.

"""
function is_equality(constraint::JuMP.ConstraintRef)::Bool
    model = constraint.model
    index = constraint.index
    set = _get_set_of_constraint(model, constraint, index)
    return set_implies_equality(set)
end

"""
    is_inequality(constraint::JuMP.ConstraintRef)

Detect whether a constraint is an inequality constraint.
"""
function is_inequality(constraint::JuMP.ConstraintRef)::Bool
    model = constraint.model
    index = constraint.index
    set = _get_set_of_constraint(model, constraint, index)
    return set_implies_inequality(set)
end

function get_equality_constraints(
    constraints::Vector{<:JuMP.ConstraintRef}
)::Vector{JuMP.ConstraintRef}
    eq_cons = Vector{JuMP.ConstraintRef}()
    for con in constraints
        if is_equality(con)
            push!(eq_cons, con)
        end
    end
    return eq_cons
end

"""
    get_equality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}

Return a vector of equality constraints in the provided model.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, v);

julia> @constraint(m, v == 1);

julia> eq_cons = MathProgIncidence.get_equality_constraints(m);

julia> display(eq_cons)
1-element Vector{ConstraintRef}:
 eq_con_1 : v = 1.0

```

"""
function get_equality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}
    constraints = JuMP.all_constraints(
        model,
        include_variable_in_set_constraints=true,
    )
    return get_equality_constraints(constraints)
end

function get_inequality_constraints(
    constraints::Vector{<:JuMP.ConstraintRef}
)::Vector{JuMP.ConstraintRef}
    ineq_cons = Vector{JuMP.ConstraintRef}()
    for con in constraints
        if is_inequality(con)
            push!(ineq_cons, con)
        end
    end
    return ineq_cons
end

"""
    get_inequality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}

Return the inequality constraints in the provided model.

# Example
```julia-repl
julia> using JuMP

julia> import MathProgIncidence

julia> m = Model();

julia> @variable(m, x[1:2] >= 0);

julia> @constraint(m, x[1]*x[2] == 1);

julia> @constraint(m, x[1] + x[2] >= 4);

julia> MathProgIncidence.get_inequality_constraints(m)
3-element Vector{ConstraintRef}:
 x[1] + x[2] ≥ 4
 x[1] ≥ 0
 x[2] ≥ 0
```
Note that variable-in-set constraints *are* included.
"""
function get_inequality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}
    constraints = JuMP.all_constraints(
        model,
        include_variable_in_set_constraints = true,
    )
    return get_inequality_constraints(constraints)
end

function get_active_inequality_constraints(
    model::JuMP.Model;
    tolerance::Float64=0.0,
)::Vector{JuMP.ConstraintRef}
    active_ineq = Vector{JuMP.ConstraintRef}()
    constraints = JuMP.all_constraints(
        model,
        include_variable_in_set_constraints = true,
    )
    for con in get_inequality_constraints(constraints)
        if is_active(con; tolerance = tolerance)
            push!(active_ineq, con)
        end
    end
    return active_ineq
end

"""
    is_active(con::JuMP.ConstraintRef; tolerance = 0.0)

Return whether the constraint is active within tolerance.

# Methods
`is_active` is supported for constraints with the following set types:
- `MOI.GreaterThan`
- `MOI.LessThan`
- `MOI.Interval`

"""
function is_active(
    con::JuMP.ConstraintRef;
    tolerance::Float64=0.0,
)::Bool
    model = con.model
    index = con.index
    set = _get_set_of_constraint(model, con, index)
    return is_active(con, set; tolerance = tolerance)
end

function is_active(
    con::JuMP.ConstraintRef,
    set::MOI.GreaterThan;
    tolerance::Float64=0.0,
)::Bool
    return abs(JuMP.value(con) - set.lower) <= tolerance
end

function is_active(
    con::JuMP.ConstraintRef,
    set::MOI.LessThan;
    tolerance::Float64=0.0,
)::Bool
    return abs(set.upper - JuMP.value(con)) <= tolerance
end

function is_active(
    con::JuMP.ConstraintRef,
    set::MOI.Interval;
    tolerance::Float64=0.0
)::Bool
    return (
        abs(set.upper - JuMP.value(con)) <= tolerance
        || abs(JuMP.value(con) - set.lower) <= tolerance
    )
end

function is_active(
    con::JuMP.ConstraintRef,
    set::MOI.AbstractSet;
    tolerance::Float64=0.0
)
    # Note that this is not a TypeError as we don't have a single that is_active
    # supports.
    throw(ArgumentError("is_active is only supported for inequality constraints"))
end

function _get_constraints(
    model::JuMP.Model;
    include_inequality::Bool,
    include_active_inequalities::Bool,
    tolerance::Float64 = 0.0,
)::Vector{JuMP.ConstraintRef}
    if include_inequality && include_active_inequalities
        throw(ArgumentError(
            "include_inequality and include_active_inequalities cannot both be true"
        ))
    end
    eq_constraints = get_equality_constraints(model)
    if include_inequality
        ineq_constraints = get_inequality_constraints(model)
        constraints = cat(eq_constraints, ineq_constraints, dims = 1)
    elseif include_active_inequalities
        ineq_constraints = get_active_inequality_constraints(
            model, tolerance = tolerance
        )
        constraints = cat(eq_constraints, ineq_constraints, dims = 1)
    else
        constraints = eq_constraints
    end
    return constraints
end

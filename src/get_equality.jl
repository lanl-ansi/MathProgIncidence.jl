#  ___________________________________________________________________________
#
#  JuMPIn.jl: JuMP Incidence Graph Analysis
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

function set_implies_equality(set::T)::Bool where T<:MOI.AbstractVectorSet
    throw(TypeError(set, MOI.AbstractScalarSet, typeof(set)))
end

# NOTE: Have not implemented tolerance arg in any calling function yet.
function set_implies_equality(
    set::MOI.Interval,
    tolerance::Float64=Float64(0)
)::Bool
    return abs(set.upper - set.lower) <= tolerance
end

# NOTE: Overload this function to define other constraint sets (e.g. Zeros)
# as equalities.
#
# I intend this to catch anything that is a subtype of ConstraintSet.
# Is the "subtype syntax" necessary?
# It appears, e.g., MOI.GreaterThan is not a subtype of MOI.ConstraintSet.
# This is non-intuitive for me...
# So this won't catch arbitrary MOI Sets...
# The right type to use is actually AbstractSet?
# Nope. This doesn't catch GreatherThan/LessThan
# ^ Wrong. Actually it does, but my syntax was wrong. I guess providing
# Type{<:SomeType} enforces that the argument must be a Type, rather than
# an instance of the type. "issubclass", where what I'm doing below is
# "isinstance"
"""
    set_implies_equality(set::T)::Bool where T<:MathOptInterface.AbstractSet

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
function set_implies_equality(
    set::T
)::Bool where T<:MOI.AbstractSet
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

function set_implies_inequality(set::MOI.AbstractSet)::Bool
    return false
end

# TODO: set_implies_inequality(set::MOI.AbstractSet, tolerance)
# that just calls set_implies_inequality(set). I.e. the default
# is to ignore the tolerance.

"""
    is_equality(constraint::JuMP.Model)::Bool

Detect whether a constraint is an equality constraint.

"""
function is_equality(constraint::JuMP.ConstraintRef)::Bool
    model = constraint.model
    index = constraint.index
    set = _get_set_of_constraint(model, constraint, index)
    return set_implies_equality(set)
end

"""
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
julia> import JuMPIn as ji
julia> m = Model();
julia> @variable(m, v);
julia> @constraint(m, v == 1);
julia> eq_cons = ji.get_equality_constraints(m);
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

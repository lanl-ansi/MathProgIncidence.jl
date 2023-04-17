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
)::Bool where T<:moi.AbstractSet
    return false
end

"""
    is_equality(constraint::JuMP.Model)::Bool

Detect whether a constraint is an equality constraint.

"""
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

"""
    get_equality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}

Return a vector of equality constraints in the provided model.

This function is also accessible via the `JuMPIn` module.

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
function get_equality_constraints(model::jmp.Model)::Vector{jmp.ConstraintRef}
    constraints = jmp.all_constraints(
        model,
        include_variable_in_set_constraints=false,
    )
    return get_equality_constraints(constraints)
end

end # module get_equality

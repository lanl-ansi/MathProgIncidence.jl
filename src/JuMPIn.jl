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

module JuMPIn
# TODO: Explicitly export functions?

# Methods to identify variables
include("identify_variables.jl") # IdentifyVariables
identify_unique_variables = IdentifyVariables.identify_unique_variables

# Methods to identify equality constraints
include("get_equality.jl") # GetEquality
get_equality_constraints = GetEquality.get_equality_constraints

# Methods to get the incidence graph of constraints and variables
include("incidence_graph.jl") # IncidenceGraph
get_bipartite_incidence_graph = IncidenceGraph.get_bipartite_incidence_graph

# Methods to apply graph algorithms to JuMP models
include("interface.jl") # Interface
IncidenceGraphInterface = Interface.IncidenceGraphInterface
get_adjacent = Interface.get_adjacent
maximum_matching = Interface.maximum_matching
dulmage_mendelsohn = Interface.dulmage_mendelsohn

end

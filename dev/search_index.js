var documenterSearchIndex = {"docs":
[{"location":"reference/index.html#Reference-Outline","page":"Reference Outline","title":"Reference Outline","text":"","category":"section"},{"location":"reference/index.html","page":"Reference Outline","title":"Reference Outline","text":"Pages = [\n    \"get_equality.md\",\n    \"identify_variables.md\",\n    \"incidence_graph.md\",\n    \"interface.md\",\n    \"incidence_matrix.md\",\n]","category":"page"},{"location":"reference/get_equality.html#Identifying-equality-and-inequality-constraints","page":"Identifying equality and inequality constraints","title":"Identifying equality and inequality constraints","text":"","category":"section"},{"location":"reference/get_equality.html","page":"Identifying equality and inequality constraints","title":"Identifying equality and inequality constraints","text":"CurrentModule = MathProgIncidence","category":"page"},{"location":"reference/get_equality.html","page":"Identifying equality and inequality constraints","title":"Identifying equality and inequality constraints","text":"get_equality_constraints\nget_inequality_constraints\nis_equality\nis_inequality\nset_implies_equality\nset_implies_inequality\nis_active","category":"page"},{"location":"reference/get_equality.html#MathProgIncidence.get_equality_constraints","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.get_equality_constraints","text":"get_equality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}\n\nReturn a vector of equality constraints in the provided model.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v);\n\njulia> @constraint(m, v == 1);\n\njulia> eq_cons = MathProgIncidence.get_equality_constraints(m);\n\njulia> display(eq_cons)\n1-element Vector{ConstraintRef}:\n eq_con_1 : v = 1.0\n\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.get_inequality_constraints","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.get_inequality_constraints","text":"get_inequality_constraints(model::JuMP.Model)::Vector{JuMP.ConstraintRef}\n\nReturn the inequality constraints in the provided model.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, x[1:2] >= 0);\n\njulia> @constraint(m, x[1]*x[2] == 1);\n\njulia> @constraint(m, x[1] + x[2] >= 4);\n\njulia> MathProgIncidence.get_inequality_constraints(m)\n3-element Vector{ConstraintRef}:\n x[1] + x[2] ≥ 4\n x[1] ≥ 0\n x[2] ≥ 0\n\nNote that variable-in-set constraints are included.\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.is_equality","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.is_equality","text":"is_equality(constraint::JuMP.ConstraintRef)::Bool\n\nDetect whether a constraint is an equality constraint.\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.is_inequality","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.is_inequality","text":"is_inequality(constraint::JuMP.ConstraintRef)\n\nDetect whether a constraint is an inequality constraint.\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.set_implies_equality","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.set_implies_equality","text":"set_implies_equality(set::MOI.AbstractSet)\n\nDetect whether the set defines an equality constraint, i.e. is a singleton.\n\nImplementation\n\nMethods are defined for the following MathOptInterface.Sets:\n\nMathOptInterface.EqualTo \nMathOptInterface.Interval\n\nIf a MathOptInterface.AbstractVectorSet is provided, an error is raised. For any other type of set, false is returned. To support additional types of constraints in is_equality and get_equality_constraints, additional methods of set_implies_equality should be defined.\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.set_implies_inequality","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.set_implies_inequality","text":"set_implies_inequality(set::MOI.AbstractSet)\n\nDetect whether the set defines an inequality constraint.\n\nThis function is defined for scalar sets. Calling with a vector set will result in a TypeError.\n\n\n\n\n\n","category":"function"},{"location":"reference/get_equality.html#MathProgIncidence.is_active","page":"Identifying equality and inequality constraints","title":"MathProgIncidence.is_active","text":"is_active(con::JuMP.ConstraintRef; tolerance = 0.0)\n\nReturn whether the constraint is active within tolerance.\n\nMethods\n\nis_active is supported for constraints with the following set types:\n\nMOI.GreaterThan\nMOI.LessThan\nMOI.Interval\n\n\n\n\n\n","category":"function"},{"location":"reference/incidence_matrix.html#Incidence-matrix","page":"Incidence matrix","title":"Incidence matrix","text":"","category":"section"},{"location":"reference/incidence_matrix.html","page":"Incidence matrix","title":"Incidence matrix","text":"CurrentModule = MathProgIncidence","category":"page"},{"location":"reference/incidence_matrix.html","page":"Incidence matrix","title":"Incidence matrix","text":"incidence_matrix","category":"page"},{"location":"reference/incidence_matrix.html#MathProgIncidence.incidence_matrix","page":"Incidence matrix","title":"MathProgIncidence.incidence_matrix","text":"incidence_matrix(constraints, variables)::SparseMatrixCSC\n\nReturn a SparseMatrixCSC with entries corresponding to incidence between the provided constraints and variables.\n\nRows correspond to constraints, columns correspond to variables, and matrix entries correspond to edges in the incidence graph. All matrix entries have values of 1.0.\n\n\n\n\n\nincidence_matrix(igraph::IncidenceGraphInterface)::SparseMatrixCSC\n\nReturn the incidence matrix associated with the provided IncidenceGraphInterface.\n\n\n\n\n\n","category":"function"},{"location":"reference/interface.html#Algorithm-interfaces","page":"Algorithm interfaces","title":"Algorithm interfaces","text":"","category":"section"},{"location":"reference/interface.html","page":"Algorithm interfaces","title":"Algorithm interfaces","text":"CurrentModule = MathProgIncidence","category":"page"},{"location":"reference/interface.html","page":"Algorithm interfaces","title":"Algorithm interfaces","text":"IncidenceGraphInterface\nget_adjacent\nmaximum_matching\ndulmage_mendelsohn\nconnected_components","category":"page"},{"location":"reference/interface.html#MathProgIncidence.IncidenceGraphInterface","page":"Algorithm interfaces","title":"MathProgIncidence.IncidenceGraphInterface","text":"IncidenceGraphInterface(\n    model;\n    include_inequality = false,\n    include_active_inequalities = false,\n    tolerance = 0.0,\n)\n\nA bipartite incidence graph of JuMP constraints and variables.\n\nThis is the primary data type accepted by the algorithms implemented in the remainder of this module. This type can be instantiated with a JuMP model or a tuple of (graph, con_node_map, var_node_map), as returned by get_bipartite_incidence_graph. If a model is provided, optional arguments are the same as those provided to get_bipartite_incidence_graph.\n\nNote that the fields of this struct are private, and may change behavior in a future release without warning.\n\nExample using only equality constraints\n\nusing JuMP\nimport MathProgIncidence\nm = Model()\n@variable(m, v[1:3] >= 0)\n@constraint(m, eq_1, v[1] + v[3]^2 == 1.0)\n@NLconstraint(m, eq_2, v[1]*v[2]^1.5 == 2.0)\ngraph = MathProgIncidence.IncidenceGraphInterface(m)\n\nExample including active inequality constraints\n\nusing JuMP\nimport Ipopt\nimport MathProgIncidence\nm = Model(Ipopt.Optimizer)\n@variable(m, v[1:3] >= 0)\n@NLconstraint(m, eq_1, v[1]*v[2]^1.5 == 2.0)\n@constraint(m, ineq_1, v[1] + v[2] + v[3] >= 7)\n@objective(m, Min, v[1]^2 + 2*v[2]^2 + 3*v[3]^2)\noptimize!(m)\ngraph = MathProgIncidence.IncidenceGraphInterface(\n    m, include_active_inequalities = true, tolerance = 1e-6\n)\n\n\n\n\n\n","category":"type"},{"location":"reference/interface.html#MathProgIncidence.get_adjacent","page":"Algorithm interfaces","title":"MathProgIncidence.get_adjacent","text":"get_adjacent(\n    igraph::IncidenceGraphInterface,\n    constraint::JuMP.ConstriantRef,\n)::Vector{JuMP.VariableRef}\n\nReturn the variables adjacent to a constraint in an incidence graph.\n\n\n\n\n\nget_adjacent(\n    igraph::IncidenceGraphInterface,\n    variable::JuMP.VariableRef,\n)::Vector{JuMP.ConstraintRef}\n\nReturn the constraints adjacent to a variable in an incidence graph.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v[1:3]);\n\njulia> @constraint(m, eq_1, v[1] + v[3] == 1);\n\njulia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);\n\njulia> igraph = MathProgIncidence.IncidenceGraphInterface(m);\n\njulia> adj_cons = MathProgIncidence.get_adjacent(igraph, v[1]);\n\njulia> display(adj_cons)\n2-element Vector{ConstraintRef}:\n eq_1 : v[1] + v[3] = 1.0\n v[1] * v[2] ^ 3.0 - 2.0 = 0\n\n\n\n\n\n\n","category":"function"},{"location":"reference/interface.html#MathProgIncidence.maximum_matching","page":"Algorithm interfaces","title":"MathProgIncidence.maximum_matching","text":"maximum_matching(igraph::IncidenceGraphInterface)::Dict\n\nCompute a maximum matching of variables and constraints in the incidence graph.\n\nThe returned Dict maps JuMP ConstraintRefs to their matched VariableRefs.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v[1:3]);\n\njulia> @constraint(m, eq_1, v[1] + v[3] == 1);\n\njulia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);\n\njulia> igraph = MathProgIncidence.IncidenceGraphInterface(m);\n\njulia> matching = MathProgIncidence.maximum_matching(igraph);\n\njulia> display(matching)\nDict{ConstraintRef, VariableRef} with 2 entries:\n  v[1] * v[2] ^ 3.0 - 2.0 = 0 => v[2]\n  eq_1 : v[1] + v[3] = 1.0 => v[1]\n\n\n\n\n\n\n","category":"function"},{"location":"reference/interface.html#MathProgIncidence.dulmage_mendelsohn","page":"Algorithm interfaces","title":"MathProgIncidence.dulmage_mendelsohn","text":"dulmage_mendelsohn(igraph::IncidenceGraphInterface)\n\nReturn the Dulmage-Mendelsohn partition of variables and constraints in an incidence graph.\n\nThe returned type is a Tuple of two NamedTuples, (con_dmp, var_dmp). These NamedTuples have the following fields:\n\ncon_dmp:\n\nunderconstrained – The constraints matched with variables that can possibly be unmatched in a maximum cardinality matching\nsquare – The constraints that cannot possibly be unmatched in a maximum matching\noverconstrained – The constraints that are matched, but can possibly be unmatched in a maximum matching\nunmatched – The constraints that are not matched in the maximum matching that was found\n\nvar_dmp:\n\nunmatched – The variables that are not matched in the maximum matching that was found\nunderconstrained – The variables that can possibly be unmatched in a maximum matching\nsquare – The variables that cannot possibly be unmatched in a maximum matching\noverconstrained – The variables matched with constraints that can possibly be unmatched in a maximum cardinality matching\n\nThe Dulmage-Mendelsohn partition groups nodes in a bipartite graph into three unique subsets. In the application to constraints and variables, these may be thought of as:\n\nthe \"overconstrained subsystem\", which has more constraints than variables,\nthe \"underconstrained subsystem\", which has more variables than constraints,\nand the \"square subsystem\", which has the same number of variables as constraints\n\nIn the NamedTuples returned by this function, the constraints in the overconstrained subsystem are split into overconstrained and unmatched, while the variables in the underconstrained subsystem are split into underconstrained and unmatched. This is because it is useful to explicitly check whether there are any unmatched variables and constraints, and also useful to recover the maximum matching by zip-ing corresponding variables and constraints.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v[1:4]);\n\njulia> @constraint(m, eq_1, v[1] + v[3] == 1);\n\njulia> @NLconstraint(m, eq_2, v[1]*v[2]^3 == 2);\n\njulia> @constraint(m, eq_3, v[4]^2 == 3);\n\njulia> igraph = MathProgIncidence.IncidenceGraphInterface(m);\n\njulia> con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph);\n\njulia> # Assert that there are no unmatched constraints\n\njulia> @assert isempty(con_dmp.unmatched);\n\njulia> display(var_dmp.unmatched)\n1-element Vector{VariableRef}:\n v[3]\n\njulia> display(var_dmp.underconstrained)\n2-element Vector{VariableRef}:\n v[1]\n v[2]\n\njulia> display(con_dmp.underconstrained)\n2-element Vector{ConstraintRef}:\n eq_1 : v[1] + v[3] = 1.0\n v[1] * v[2] ^ 3.0 - 2.0 = 0\n \njulia> display(var_dmp.square)\n1-element Vector{VariableRef}:\n v[4]\n\njulia> display(con_dmp.square)\n1-element Vector{ConstraintRef}:\n eq_3 : v[4]² = 3.0\n\njulia> # As there are no unmatched constraints, the overconstrained subsystem is empty\n\n\n\n\n\n","category":"function"},{"location":"reference/interface.html#MathProgIncidence.connected_components","page":"Algorithm interfaces","title":"MathProgIncidence.connected_components","text":"connected_components(igraph::IncidenceGraphInterface)\n\nReturn the connected components of a bipartite incidence graph of constraints and variables.\n\nThe connected components are returned as two vector-of-vectors, containing the variables in each connected component and the constraints in each connected component. Note that the input graph is undirected, so there is no distinction between strongly and weakly connected components.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, x[1:2] >= 0);\n\njulia> @constraint(m, eq1, x[1] == 1);\n\njulia> @constraint(m, eq2, x[2]^2 == 2);\n\njulia> igraph = MathProgIncidence.IncidenceGraphInterface(m);\n\njulia> con_comps, var_comps = MathProgIncidence.connected_components(igraph);\n\njulia> con_comps\n2-element Vector{Vector{ConstraintRef}}:\n [eq1 : x[1] = 1]\n [eq2 : x[2]² = 2]\n\njulia> var_comps\n2-element Vector{Vector{VariableRef}}:\n [x[1]]\n [x[2]]\n\n\n\n\n\n\nconnected_components(constraints, variables)\n\nReturn the connected components of a bipartite incidence graph of constraints and variables.\n\nThe method that accepts constraints and variables directly is convenient for working with the output of the Dulmage-Mendelsohn partition. It is often used to decompose and help debug the over and under-constrained subsystems.\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, x[1:4] >= 0);\n\njulia> @constraint(m, eq1, x[1] + x[3] == 7);\n\njulia> @constraint(m, eq2, x[2]^2 + x[4]^2 == 1);\n\njulia> igraph = MathProgIncidence.IncidenceGraphInterface(m);\n\njulia> con_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph);\n\njulia> uc_con = con_dmp.underconstrained;\n\njulia> uc_var = [var_dmp.unmatched..., var_dmp.underconstrained...];\n\njulia> con_comps, var_comps = MathProgIncidence.connected_components(uc_con, uc_var);\n\njulia> con_comps\n2-element Vector{Vector{ConstraintRef}}:\n [eq1 : x[1] + x[3] = 7]\n [eq2 : x[2]² + x[4]² = 1]\n\njulia> var_comps\n2-element Vector{Vector{VariableRef}}:\n [x[3], x[1]]\n [x[4], x[2]]\n\n\n\n\n\n\n","category":"function"},{"location":"overview.html#Overview","page":"Overview","title":"Overview","text":"","category":"section"},{"location":"overview.html#What-is-MathProgIncidence?","page":"Overview","title":"What is MathProgIncidence?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"MathProgIncidence is a JuMP extension that provides algorithms for analyzing incidence graphs or matrices defined by JuMP variables and constraints. The standard graph that is analyzed is a bipartite graph of variables and constraints, where an edge exists between a variable and constraint if the variable participates in the constraint.","category":"page"},{"location":"overview.html#Why-is-MathProgIncidence-useful?","page":"Overview","title":"Why is MathProgIncidence useful?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"In a large modeling/optimization project, especially one involving several developers, it is fairly easy to make a mistake designing or implementing an algebraic model. These mistakes commonly cause singularities in the Jacobian of equality constraints. The algorithms implemented in MathProgIncidence allow a modeler to identify irreducible subsets of variables and constraints that are causing singularities, which can be very useful when debugging a suspected modeling error.","category":"page"},{"location":"overview.html#When-should-I-suspect-a-modeling-error?","page":"Overview","title":"When should I suspect a modeling error?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"One obvious situation is that the solution of optimization problems yields solutions that do not make sense. Another is that optimization problems fail to converge in a reasonable amount of time. However, optimization problems can fail to converge for a wide variety of reasons, and the reason for non-convergence is highly dependent on the algorithm used and optimization problem being solved. In the case of nonlinear local optimization, symptoms that are often indicative of modeling errors are large regularization coefficients and large numbers of restoration iterations in interior point methods.","category":"page"},{"location":"overview.html#What-algorithms-does-MathProgIncidence-implement?","page":"Overview","title":"What algorithms does MathProgIncidence implement?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"The Dulmage-Mendelsohn partition, which detects subsets of variables and constraints causing a structural singularity. See MathProgIncidence.dulmage_mendelsohn.","category":"page"},{"location":"overview.html","page":"Overview","title":"Overview","text":"More algorithms may be implemented in the future.","category":"page"},{"location":"overview.html#What-models-should-these-be-applied-to?","page":"Overview","title":"What models should these be applied to?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"These algorithms should be used for applications where singularity of (a subsystem of) constraints and variables violates an assumption made by the mathematical/physical model or the solver. Examples include:","category":"page"},{"location":"overview.html","page":"Overview","title":"Overview","text":"Index-1 differential-algebraic equations (DAEs)\nChemical process flowsheets\nGas pipeline networks","category":"page"},{"location":"overview.html","page":"Overview","title":"Overview","text":"Typically, models of highly nonlinear phenomena distributed over some network (or time period) are error-prone to implement and could benefit from these algorithms. In addition, nonlinear local optimization algorithms such as Ipopt (which are used as subroutines of global and mixed-integer nonlinear solvers) assume that the Jacobian of equality constraints is full row rank. A convenient way to check this is to fix the variables that you intend to be degrees of freedom, then check the Jacobian of equality constraints for singularity. If this Jacobian is nonsingular, the full equality Jacobian is full row rank. Otherwise, this assumption may be violated, and these algorithms may help deteremine the reason why.","category":"page"},{"location":"overview.html#How-can-I-cite-MathProgIncidence?","page":"Overview","title":"How can I cite MathProgIncidence?","text":"","category":"section"},{"location":"overview.html","page":"Overview","title":"Overview","text":"If you use MathProgIncidence in your research, we would appreciate you citing the following paper:","category":"page"},{"location":"overview.html","page":"Overview","title":"Overview","text":"@article{parker2023dulmage,\ntitle = {Applications of the {Dulmage-Mendelsohn} decomposition for debugging nonlinear optimization problems},\njournal = {Computers \\& Chemical Engineering},\nvolume = {178},\npages = {108383},\nyear = {2023},\nissn = {0098-1354},\ndoi = {https://doi.org/10.1016/j.compchemeng.2023.108383},\nurl = {https://www.sciencedirect.com/science/article/pii/S0098135423002533},\nauthor = {Robert B. Parker and Bethany L. Nicholson and John D. Siirola and Lorenz T. Biegler},\n}","category":"page"},{"location":"reference/identify_variables.html#Identifying-variables-in-constraints","page":"Identifying variables in constraints","title":"Identifying variables in constraints","text":"","category":"section"},{"location":"reference/identify_variables.html","page":"Identifying variables in constraints","title":"Identifying variables in constraints","text":"CurrentModule = MathProgIncidence","category":"page"},{"location":"reference/identify_variables.html","page":"Identifying variables in constraints","title":"Identifying variables in constraints","text":"identify_unique_variables\n_get_variable_terms","category":"page"},{"location":"reference/identify_variables.html#MathProgIncidence.identify_unique_variables","page":"Identifying variables in constraints","title":"MathProgIncidence.identify_unique_variables","text":"identify_unique_variables(constraints::Vector)::Vector{JuMP.VariableRef}\n\nReturn a vector of variables that participate in the provided constraints.\n\nEach variable appears at most one time in the returned vector. If we receive a vector, we just assume it is a vector of constraints. This is because I couldn't get an argument of type Vector{ConstraintRef} to work...\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v[1:3]);\n\njulia> @constraint(m, eq_1, v[2] == 1);\n\njulia> @NLconstraint(m, eq_2, v[2]*v[3]^1.5 == 2);\n\njulia> vars = MathProgIncidence.identify_unique_variables([eq_1, eq_2]);\n\njulia> display(vars)\n2-element Vector{VariableRef}:\n v[2]\n v[3]\n\n\n\n\n\n\nidentify_unique_variables(\n    model::JuMP.Model,\n    include_inequality = false,\n    include_active_inequalities = false,\n    tolerance = 0.0,\n)\n\nReturn a vector of variables that participate in constraints in the model.\n\nEach variable appears at most one time in the returned vector.\n\nOptional keyword arguments\n\ninclude_inequality: Whether to include variables that participate in any inequality constraints\ninclude_active_inequalities: Whether to include variables that participate active inequality constraints (at the latest solution). include_active_inequalities and include_inequality are mutually exclusive.\ntolerance: Tolerance used to determine if an inequality constraint is active\n\n\n\n\n\nidentify_unique_variables(constraint::JuMP.ConstraintRef)\n\nReturn a vector containing the variables that participate in this constraint.\n\nEach variable appears at most one time in the returned vector.\n\n\n\n\n\nidentify_unique_variables(\n    constraint::JuMP.ConstraintRef,\n    index::JuMP.ConstraintIndex,\n)::Vector{JuMP.VariableRef}\n\nReturn a vector containing the variables that participate in this constraint.\n\nEach variable appears at most one time in the returned vector. The index argument is provided so we know whether we have a \"regular\" constraint or a nonlinear constraint.\n\n\n\n\n\nidentify_unique_variables(\n    constraint::JuMP.ConstraintRef,\n    index::MathOptInterface.Nonlinear.ConstraintIndex,\n)::Vector{JuMP.VariableRef}\n\nReturn a vector containing the variables that participate in this constraint.\n\nEach variable appears at most one time in the returned vector. The index argument is provided so we know whether we have a \"regular\" constraint or a nonlinear constraint.\n\n\n\n\n\nidentify_unique_variables(fcn)::Vector{JuMP.VariableIndex}\n\nReturn the variables that appear in the provided MathOptInterface function.\n\nNo variable will appear more than once.\n\nImplementation\n\nOnly ScalarNonlinearFunction, ScalarQuadraticFunction, and ScalarAffineFunction are supported. This can be changed there is demand for other functions. These methods are implemented by first identifying all variables that participate in the function, then filtering out duplicate variables.\n\n\n\n\n\nidentify_unique_variables(term)\n\nReturn the variables that participate in the provided term of a MathOptInterface function.\n\nA variable appears at most one time in the returned vector.\n\nImplementation\n\nCurrently implemented only for ScalarAffineTerm and ScalarQuadraticTerm.\n\n\n\n\n\n","category":"function"},{"location":"reference/identify_variables.html#MathProgIncidence._get_variable_terms","page":"Identifying variables in constraints","title":"MathProgIncidence._get_variable_terms","text":"_get_variable_terms(fcn)\n\nReturn a tuple of vectors of terms for the provided MathOptInterface function.\n\nImplementation\n\nCurrently implemented only for ScalarQuadraticFunction and ScalarAffineFunction.\n\n\n\n\n\n","category":"function"},{"location":"reference/incidence_graph.html#Incidence-graph","page":"Incidence graph","title":"Incidence graph","text":"","category":"section"},{"location":"reference/incidence_graph.html","page":"Incidence graph","title":"Incidence graph","text":"CurrentModule = MathProgIncidence","category":"page"},{"location":"reference/incidence_graph.html","page":"Incidence graph","title":"Incidence graph","text":"get_bipartite_incidence_graph","category":"page"},{"location":"reference/incidence_graph.html#MathProgIncidence.get_bipartite_incidence_graph","page":"Incidence graph","title":"MathProgIncidence.get_bipartite_incidence_graph","text":"get_bipartite_incidence_graph(\n    model,\n    include_inequality = false,\n    include_active_inequalities = false,\n    tolerance = 0.0,\n)\n\nReturn the bipartite incidence graph of (scalar) variables and constraints in the JuMP model.\n\nArguments\n\nmodel::JuMP.Model: Model whose incidence graph is constructed\ninclude_inequality::Bool: Whether all inequality constraints should be included in the graph\ninclude_active_inequalities::Bool: Whether only active inequality constraints should be included in the graph. Constraint activity is determined using the most recent solution. include_active_inequalities and include_inequality are mutually exclusive.\ntolerance::Float64: Tolerance for deciding whether an inequality constraint is active\n\nReturns\n\nThis function returns a tuple (graph, con_node_map, var_node_map)\n\ngraph – a tuple (A, B, E) where A and B contain the integer nodes in the bipartite sets for constraints and variables and E contains edges in the form of tuples of integers (a, b), where a is in  A and b is in B.\ncon_node_map – a Dict mapping JuMP ConstraintRefs to nodes\nvar_node_map – a Dict mapping JuMP VariableRefs to nodes\n\nThe constraints in the graph are all the (by default, equality) constraints in the model, and the variables are those that participate in these constraints (not all variables in the model).\n\nExample\n\njulia> using JuMP\n\njulia> import MathProgIncidence\n\njulia> m = Model();\n\njulia> @variable(m, v[1:3]);\n\njulia> @constraint(m, eq_1, v[1] + v[3]^2 == 1.0);\n\njulia> @NLconstraint(m, eq_2, v[1]*v[2]^1.5 == 2.0);\n\njulia> graph, con_node_map, var_node_map = MathProgIncidence.get_bipartite_incidence_graph(m);\n\njulia> A, B, E = graph;\n\njulia> M = length(A);\n\njulia> N = length(B);\n\njulia> imat = zeros(M, N);\n\njulia> for (a, b) in E\n           imat[a, b-M] = 1.0;\n       end\n\njulia> display(imat)\n2×3 Matrix{Float64}:\n 1.0  1.0  0.0\n 0.0  1.0  1.0\n\n\nConvention\n\nThe returned graph follows a convention where, for a JuMP model with N constraints and M variables, the first N nodes are constraints and nodes N+1 through N+M are variables. Nodes are always contiguous integers starting at 1.\n\nMethods\n\nMethods are implemented that accept a JuMP Model, a vector of ConstraintRefs, and vectors of ConstraintRefs and VariableRefs. If variables are provided, then only these variables participate in the graph, regardless of which variables participate in the constraints.\n\n\n\n\n\n","category":"function"},{"location":"index.html#MathProgIncidence.jl","page":"Introduction","title":"MathProgIncidence.jl","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Math Programming Incidence graph analysis. Tools for constructing and analyzing the incidence graph or matrix of variables and constraints in a JuMP model.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"If you are wondering what JuMP is, please see here. If you are familiar with JuMP and wondering whether you should use MathProgIncidence, you may find the Overview helpful.","category":"page"},{"location":"index.html#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"MathProgIncidence is registered on the Julia general package registry and can be installed with:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"julia> ]\n(@v1.X) pkg> add MathProgIncidence","category":"page"},{"location":"index.html#Resources","page":"Introduction","title":"Resources","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"This documentation contains the following resources:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Simple Example – A simple demonstration of the intended use of MathProgIncidence\nReference Outline – Call signatures and documentation for the public API of MathProgIncidence","category":"page"},{"location":"example.html#Simple-Example","page":"Simple Example","title":"Simple Example","text":"","category":"section"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"This page walks through a simple example of using the Dulmage-Mendelsohn partition to debug a structural singularity.","category":"page"},{"location":"example.html#Dulmage-Mendelsohn","page":"Simple Example","title":"Dulmage-Mendelsohn","text":"","category":"section"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"We start with some imports and by creating a JuMP model we are interested in. Usually the model we are interested in debugging is much larger and more complicated than this. This particular system appeared when debugging a dynamic 1-D partial differential-algebraic equation (PDAE) model representing a chemical looping combustion reactor.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"using JuMP\nimport MathProgIncidence\n\nm = Model()\ncomps = [1, 2, 3]\n@variable(m, x[comps], start=1/3.0)\n@variable(m, flow_comp[comps], start=10.0)\n@variable(m, flow, start=30.0)\n@variable(m, rho, start=1.0)\n\n@constraint(m, sum_comp_eqn, sum(x) == 1)\n@constraint(m, comp_dens_eqn, x*rho .== [1.0, 1.1, 1.2])\n@NLconstraint(m, bulk_dens_eqn, 1/rho - sum(1/x[j] for j in comps) == 0)\n@constraint(m, comp_flow_eqn, x.*flow .== flow_comp)","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"To check this model for structural singularity, we apply the Dulmage-Mendelsohn partition.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"igraph = MathProgIncidence.IncidenceGraphInterface(m)\ncon_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"If any variables or constraints are unmatched, the (Jacobian of the) model is structurally singular.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"println(\"Unmatched constraints\")\nprintln(\"---------------------\")\nfor con in con_dmp.unmatched\n    println(\"  $con\")\nend\nprintln(\"Unmatched variables\")\nprintln(\"-------------------\")\nfor var in var_dmp.unmatched\n    println(\"  $var\")\nend","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"Unmatched constraints\n---------------------\n(1.0 / rho - (1.0 / x[1] + 1.0 / x[2] + 1.0 / x[3])) - 0.0 = 0\nUnmatched variables\n-------------------\nflow_comp[1]","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"This model has one unmatched constraint and one unmatched variable, so it is structurally singular. However, the unmatched constraint and variable are not unique. For example, flow_comp[2] could have been unmatched instead of flow_comp[1].","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"Unique subsets of variables and constraints that are useful when debugging a structural singularity are the underconstrained and overconstrained subsystems. The variables in the underconstrained subsystems are contained in the unmatched and underconstrained fields, while the constraints are contained in the underconstrained field. The variables in the overconstrained subsystem are contained in the overconstrained field, while the constraints are contained in the overconstrained and unmatched fields.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"We now construct the underconstrained and overconstrained subsystems.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"oc_con = cat(con_dmp.overconstrained, con_dmp.unmatched, dims = 1)\noc_var = var_dmp.overconstrained\nuc_con = con_dmp.underconstrained\nuc_var = cat(var_dmp.unmatched, var_dmp.underconstrained, dims = 1)","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"And display the constraints and variables contained in each.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"println(\"Overconstrained subsystem\")\nprintln(\"-------------------------\")\nprintln(\"Constraints\")\nfor con in oc_con\n    println(\"  $con\")\nend\nprintln(\"Variables\")\nfor var in oc_var\n    println(\"  $var\")\nend\nprintln()\n\nprintln(\"Underconstrained subsystem\")\nprintln(\"--------------------------\")\nprintln(\"Constraints\")\nfor con in uc_con\n    println(\"  $con\")\nend\nprintln(\"Variables\")\nfor var in uc_var\n    println(\"  $var\")\nend","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"Overconstrained subsystem\n--------------------------\nConstraints\n  sum_comp_eqn : x[1] + x[2] + x[3] = 1.0\n  comp_dens_eqn : x[1]*rho = 1.0\n  comp_dens_eqn : x[2]*rho = 1.1\n  comp_dens_eqn : x[3]*rho = 1.2\n  (1.0 / rho - (1.0 / x[1] + 1.0 / x[2] + 1.0 / x[3])) - 0.0 = 0\nVariables\n  x[1]\n  rho\n  x[2]\n  x[3]\n\nUnderconstrained subsystem\n--------------------------\nConstraints\n  comp_flow_eqn : x[1]*flow - flow_comp[1] = 0.0\n  comp_flow_eqn : x[2]*flow - flow_comp[2] = 0.0\n  comp_flow_eqn : x[3]*flow - flow_comp[3] = 0.0\nVariables\n  flow_comp[1]\n  flow\n  flow_comp[2]\n  flow_comp[3]","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"At this point we must use our intuition about the system being modeled to identify \"what is causing\" the singularity. Looking at the under and over-constrained systems, it appears that we are missing an equation to calculate flow, the total flow rate, and that rho, the bulk density, is over-specified as it is computed by both the bulk density equation and one of the component density equations.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"With this knowledge, we can eventually figure out that (a) we need an equation to calculate flow from density, and that (b) our \"bulk density equation\" is actually a skeletal density equation. Admittedly, this is difficult to figure out without the full context behind this particular system.","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"The following script constructs a new version of the model and checks it for structural singularity:","category":"page"},{"location":"example.html","page":"Simple Example","title":"Simple Example","text":"using JuMP\nimport MathProgIncidence\n\nm = Model()\ncomps = [1, 2, 3]\n@variable(m, x[comps], start=1/3.0)\n@variable(m, flow_comp[comps], start=10.0)\n@variable(m, flow, start=30.0)\n@variable(m, rho_bulk, start=1.0)\n@variable(m, rho_skel, start=1.0)\n@variable(m, porosity, start=0.25)\nvelocity = 1.0\n\n@constraint(m, sum_comp_eqn, sum(x) == 1)\n@constraint(m, comp_dens_eqn, x*rho_bulk .== [1.0, 1.1, 1.2])\n@NLconstraint(m, skel_dens_eqn, 1/rho_skel - sum(1/x[j] for j in comps) == 0)\n@constraint(m, bulk_dens_eqn, rho_bulk == (1 - porosity)*rho_skel)\n@constraint(m, flow_eqn, flow == velocity * rho_bulk)\n@constraint(m, comp_flow_eqn, x.*flow .== flow_comp)\n\nigraph = MathProgIncidence.IncidenceGraphInterface(m)\ncon_dmp, var_dmp = MathProgIncidence.dulmage_mendelsohn(igraph)\n@assert isempty(con_dmp.unmatched) && isempty(var_dmp.unmatched)","category":"page"}]
}

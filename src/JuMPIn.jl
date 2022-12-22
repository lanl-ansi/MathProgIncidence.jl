module JuMPIn
# TODO: Explicitly export functions?

# Methods to identify variables
include("identify_variables.jl") # IdentifyVariables
identify_unique_variables = IdentifyVariables.identify_unique_variables

end

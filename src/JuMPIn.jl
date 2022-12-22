module JuMPIn

# Write your package code here.

# Methods to identify variables
include("identify_variables.jl") # IdentifyVariables
identify_variables = IdentifyVariables.identify_unique_variables

end

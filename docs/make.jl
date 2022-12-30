using Documenter
using JuMPIn

makedocs(
    sitename = "JuMPIn",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "Overview" => "overview.md",
        "Simple Example" => "example.md",
        "API Reference" => [
            "reference/get_equality.md",
            "reference/identify_variables.md",
            "reference/incidence_graph.md",
            "reference/interface.md",
        ],
    ],
)

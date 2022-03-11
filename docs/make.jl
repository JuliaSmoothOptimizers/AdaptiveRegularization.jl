ENV["GKSwstype"] = "100"
using ADNLPModels
using Documenter
using Printf
using ARCTR

pages = [
    "Introduction" => "index.md",
    "Tutorial" => "benchmark.md",
    "Reference" => "reference.md",
]

makedocs(
    sitename = "ARCTR.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [ARCTR],
    pages = pages,
)

deploydocs(repo = "github.com/tmigot/ARCTR.jl.git", push_preview = true, devbranch = "main")

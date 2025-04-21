using AdaptiveRegularization
using Documenter

DocMeta.setdocmeta!(AdaptiveRegularization, :DocTestSetup, :(using AdaptiveRegularization); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [AdaptiveRegularization],
    authors = "Dussault, Jean-Pierre <Jean-Pierre.Dussault@usherbrooke.ca >, Goyette, Samuel <samuel.goyette@usherbrooke.ca>, Migot, Tangi <tangi.migot@gmail.com>",
    repo = "https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/blob/{commit}{path}#{line}",
    sitename = "AdaptiveRegularization.jl",
    format = Documenter.HTML(; canonical = "https://JuliaSmoothOptimizers.github.io/AdaptiveRegularization.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl")

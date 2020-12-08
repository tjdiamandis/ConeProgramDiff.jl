using ConeProgramDiff
using Documenter

makedocs(;
    modules=[ConeProgramDiff],
    authors="Chandler Squires, Theo Diamandis",
    repo="https://github.com/tjdiamandis/ConeProgramDiff.jl/blob/{commit}{path}#L{line}",
    sitename="ConeProgramDiff.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tjdiamandis.github.io/ConeProgramDiff.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tjdiamandis/ConeProgramDiff.jl",
)

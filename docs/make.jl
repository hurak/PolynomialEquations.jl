using PolynomialEquations
using Documenter

makedocs(;
    modules=[PolynomialEquations],
    authors="Zdenek Hurak",
    repo="https://github.com/hurak/PolynomialEquations.jl/blob/{commit}{path}#L{line}",
    sitename="PolynomialEquations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hurak.github.io/PolynomialEquations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hurak/PolynomialEquations.jl",
)

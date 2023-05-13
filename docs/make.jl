using LinearNoiseApproximation
using Documenter

DocMeta.setdocmeta!(LinearNoiseApproximation, :DocTestSetup, :(using LinearNoiseApproximation); recursive=true)

makedocs(;
    modules=[LinearNoiseApproximation],
    authors="Xiaoming Fu",
    repo="https://github.com/palmtree2013/LinearNoiseApproximation.jl/blob/{commit}{path}#{line}",
    sitename="LinearNoiseApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palmtree2013.github.io/LinearNoiseApproximation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/palmtree2013/LinearNoiseApproximation.jl",
    devbranch="main",
)
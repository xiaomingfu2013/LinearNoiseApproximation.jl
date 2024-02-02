using LinearNoiseApproximation
using Documenter

DocMeta.setdocmeta!(
    LinearNoiseApproximation,
    :DocTestSetup,
    :(using LinearNoiseApproximation);
    recursive=true,
)

makedocs(;
    modules=[LinearNoiseApproximation],
    authors="Xiaoming Fu",
    repo="https://github.com/xiaomingfu2013/LinearNoiseApproximation.jl/blob/{commit}{path}#{line}",
    sitename="LinearNoiseApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xiaomingfu2013.github.io/LinearNoiseApproximation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md", "Tutorial" => "PredatorPreyTutorial.md", "API" => "api.md"],
    warnonly=true,
)

deploydocs(; repo="github.com/xiaomingfu2013/LinearNoiseApproximation.jl", devbranch="main")

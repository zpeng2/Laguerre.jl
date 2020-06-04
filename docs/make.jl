using Laguerre
using Documenter

makedocs(;
    modules=[Laguerre],
    authors="Zhiwei Peng",
    repo="https://github.com/zpeng2/Laguerre.jl/blob/{commit}{path}#L{line}",
    sitename="Laguerre.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zpeng2.github.io/Laguerre.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zpeng2/Laguerre.jl",
)

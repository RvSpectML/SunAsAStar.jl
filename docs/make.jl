using SunAsAStar
using Documenter

DocMeta.setdocmeta!(SunAsAStar, :DocTestSetup, :(using SunAsAStar); recursive=true)

makedocs(;
    modules=[SunAsAStar],
    authors="Eric Ford, Andrea Lin, Shubham Kanodia",
    repo="https://github.com/RvSpectML/SunAsAStar.jl/blob/{commit}{path}#{line}",
    sitename="SunAsAStar.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/SunAsAStar.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RvSpectML/SunAsAStar.jl",
    devbranch="main",
)

using DFOQuadModel
using Documenter

DocMeta.setdocmeta!(DFOQuadModel, :DocTestSetup, :(using DFOQuadModel); recursive=true)

makedocs(;
    modules=[DFOQuadModel],
    authors="Iffan Hannanu <iffan.hannanu@gmail.com> and contributors",
    sitename="DFOQuadModel.jl",
    format=Documenter.HTML(;
        canonical="https://iffanh.github.io/DFOQuadModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iffanh/DFOQuadModel.jl",
    devbranch="main",
)

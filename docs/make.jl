using Documenter, EITData
using Documenter: deploydocs

makedocs(
  modules=[EITData],
  clean=true,
  pages=[
     "Home" => "index.md",
     "Tutorial" => "tutorial.md",
     "API Reference" => "api.md"
  ]
)
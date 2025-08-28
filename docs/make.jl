using Documenter

makedocs(
    sitename = "Parametric Surface Visualization - Compact Bonnet Pairs",
    authors = "Kiba Mangandango",
    pages = [
    "Home"     => "index.md",
    "Theory"   => "theory.md",
    "Examples" => "examples.md",
    "API"      => "api.md",
    ],
    # optional: configure HTML/math; leave default for now
    format = Documenter.HTML()
)

# deploy to gh-pages (replace user/repo)
deploydocs(repo = "github.com/Okami4D/Compact-Bonnet-Surface-Visualization.git")
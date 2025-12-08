# Julia Project Dependencies

This Julia project uses `Project.toml` (similar to Python's `requirements.txt`) to manage dependencies.

## Project Information

- **Name:** IEEE519
- **Version:** 0.1.0
- **Author:** Serin Thomas

## Managing Dependencies

### Adding Dependencies

To add a new package to your project:

```julia
using Pkg
Pkg.add("PackageName")
```

Or from the Julia REPL, press `]` to enter Pkg mode:

```
] add PackageName
```

### Current Dependencies

None currently specified. Add packages as needed to `Project.toml`.

## Recommended packages for this project

The `src/IEEE519.jl` script uses the following packages (some are part of Julia's standard library):

- `FFTW` (signal processing)
- `Statistics` (stdlib)
- `Plots` (plotting)
- `PlotlyJS` (optional interactive plotting backend)
- `GR` (fallback plotting backend)
- `CSV` (CSV read/write)
- `DataFrames` (tabular data)

If you want interactive Plotly plots (the script calls `plotlyjs()`), you should install `PlotlyJS` and `WebIO` and run the build step for `WebIO` so its JS runtime is installed.

## Installation (recommended)

Open the Julia REPL in the project folder and run:

```julia
using Pkg
Pkg.activate(".")        # activate this project's environment
Pkg.instantiate()         # optional: installs from Project/Manifest if present
Pkg.add(["Plots","PlotlyJS","WebIO","CSV","DataFrames","FFTW","GR"]) 
Pkg.build("WebIO")      # ensures JavaScript runtime dependencies are installed
```

Notes:

- If you don't need interactive Plotly output, you can skip `PlotlyJS` and `WebIO`.
- The script will automatically fall back to the `GR` backend if PlotlyJS fails to initialize.
- On some systems the `WebIO` build may require `nodejs`/`npm` on PATH; `Pkg.build("WebIO")` will try to install a compatible Node.js, but system policies may require manual install.

## Using the code

After installing dependencies you can run the example with:

```julia
include("src/IEEE519.jl")
```

If you see a warning about PlotlyJS / WebIO on startup, either install the missing packages above or rely on the GR backend which provides static plots.

## Dependency File Format

- **`Project.toml`** - Declares package name, version, and direct dependencies (like `setup.py`)
- **`Manifest.toml`** - Auto-generated lock file that pins exact versions of all dependencies (like `pip freeze`)

## Installation

To install all dependencies from `Project.toml`:

```julia
using Pkg
Pkg.instantiate()
```

Or in Pkg mode:

```
] instantiate
```

## Example Dependencies

When you add packages, they will appear in `Project.toml` like:

```toml
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
```

## Version Specifications

- `PackageName` - Use latest compatible version
- `PackageName@1.2` - Use version 1.2.x
- `PackageName@1.2.3` - Use exact version
- `PackageName@">=1.0,<2.0"` - Use version range

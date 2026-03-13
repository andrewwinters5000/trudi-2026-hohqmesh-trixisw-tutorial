# TRUDI 2026 -- HOHQMesh and TrixiShallowWater tutorials

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This repository contains information and code to reproduce the results presented in the
hands-on session for HOHQMesh / HOHQMesh.jl and TrixiShallowWater.jl at TRUDI 2026 in Cologne, Germany.

## Installation

To reproduce the tutorials presented herein, you need to install [Julia](https://julialang.org/).
The tutorials were prepared using Julia v1.11.3, but other versions should work, e.g., v 1.10.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface.
After cloning, you will need to instantiate and install the necessary packages.
Navigate to the main folder of this repository and execute
```shell
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```
Then, you can start a Julia REPL
```shell
julia --project=.
```
There are further instructions described in the `README.md` files in the folders
`HOHQMesh/` or `TrixiSW/` that discuss how to put together the examples.

## Authors

- Patrick Ersing (Linköping University, Sweden)
- David A. Kopriva (Florida State University, Florida, USA and San Diego State University, California, USA)
- Andrew R. Winters (Linköping University, Sweden)

## License

The code in this repository is published under the MIT license, see the
`LICENSE.md` file.

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
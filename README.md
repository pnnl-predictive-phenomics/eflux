# E-Flux (combining Expression Data with Fluxes)
[![Rye](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/rye/main/artwork/badge.json)](https://rye-up.com)

Package to run E-Flux on `cobra` models with provided transcriptomics expression data.


## Getting Started

To get started run the following:

```python
from eflux import eflux2

eflux2(cobra_model, transcriptomics)
```

## Installation

The most recent code and data can be installed directly from GitHub with:

```shell
pip install git+https://github.com/pnnl-predictive-phenomics/eflux.git
```

## License

The code in this package is licensed under the BSD-2 License.


## Development

To contribute to the eflux repo, you can install the code via:

```shell
git clone git+https://github.com/pnnl-predictive-phenomics/eflux.git
cd eflux
pip install -e .
```

Then you can create a virtual environment through:

```shell
rye sync
```

### 🥼 Testing

After cloning the repository, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
rye test
```

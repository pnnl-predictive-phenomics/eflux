# Contributing

## Issue and PRs

The first step to making a contribution is seeing if your problem is already being addressed by a current issue. If so, consider commenting to share that you are also having the same issue. 

If not, consider creating an issue to inform us of the problem, and then if you are interested in making a contribution, create a Pull Request.

## Development Environment 👷

To make a contribution, you will need to set up the development environment. This repo uses [Rye](https://github.com/astral-sh/rye), but any `virtualenv` python tool will work (Pyenv, virtualenv, venv, etc.). 

### Using Rye ✅

> Make sure to have `rye` installed before performing these steps.

If using rye, you can quickly spin up the development environment using the following:

```shell
git clone git+https://github.com/pnnl-predictive-phenomics/eflux.git
cd eflux
rye sync
```
> Note that you must first uninstall any current version you may already have.

This will will install the code in an editable configuration and install all the required dev dependencies.

### Not using Rye ❌

First you must install the code via:
```shell
git clone git+https://github.com/pnnl-predictive-phenomics/eflux.git
cd eflux
pip install -e .
```
> Note that you must first uninstall any current version you may already have.

Then you can make sure you have all the dev dependencies installed by running
```shell
pip install -r requirements-dev.lock
```

## Linting, Formatting, and Testing 🧹🎨🧪

When making commits to your branch to prepare for your PR, make sure to follow proper code quality assessments by frequently running the linter, formatter, and test suite on your code. The PR will not be able to be merged if the linter and tests do not pass. 

> Be in the root directory of the repo before running these instructions.

### Using Rye ✅

To lint your code, use the command:
```shell
rye lint
rye run mypy --ignore-missing-imports .
```
This will return a list of current issues with your code.

---
To format your code, use the command:
```shell
rye fmt
```

---
To run the test suite on your code, use the command:
```bash
rye test
```

### Not using Rye ❌

To lint your code, use the command:
```shell
ruff check .
mypy --ignore-missing-imports .
```
This will return a list of current issues with your code.

---
To format your code, use the command:
```shell
ruff format .
```

---
To test your code, use the command:
```shell
pytest
```
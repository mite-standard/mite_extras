mite_extras
==========

[![DOI](https://zenodo.org/badge/804997522.svg)](https://doi.org/10.5281/zenodo.13941745) [![PyPI version](https://badge.fury.io/py/mite-extras.svg)](https://badge.fury.io/py/mite-extras)

This package contains functionality to validate entries of the Minimum Information about a Tailoring Enzyme (MITE) repository.

This includes:

- Validation of the MITE json format (does the entry follow the `mite_schema` format)
- Validation of reaction SMARTS (does the reaction SMARTS lead to the expected product when applied to a specified substrate)
- Sanitation of SMILES and reaction SMARTS
- Fetching of NCBI GenPept/UniProtKB IDs

When used as CLI, `mite_extras` can automatically update MITE entries (see below). However, it can also be used as a library (e.g. as done in `mite_data`).

For more information, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

## Installation

**Attention: this program has only been tested with Ubuntu Linux.**

### With `pip` from PyPI

- Install with `pip install mite_extras`

### With `hatch` from GitHub

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Download or clone the [repository](https://github.com/mmzdouc/mite_extras)
- Run `hatch -v env create`

## Quick Start: Example

### Run from command line:

To validate MITE entries or update them to a new schema version (requires `mite_extras` to be installed via `pip`).

- `mite_extras -i input/ -o output/`

### Run with `hatch`:

Validate MITE entries or update them to a new schema version

- `hatch run mite_extras -i input/ -o output/`

## For devs

- Install developer dependencies with `hatch -v env create dev`
- Initialize `pre-commit` with `hatch run dev:pre-commit install`
- Run tests with `hatch run dev:pytest`
- Run CLI with `hatch run dev:mite_extras` and the appropriate options
- If necessary, remove the environment again with `hatch env remove dev`
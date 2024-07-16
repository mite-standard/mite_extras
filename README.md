mite_extras
==========

Validation and conversion functionality for reaction SMARTS collected in the scope of MITE

## Installation

### With `hatch` from GitHub
- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Download or clone the [repository](https://github.com/mmzdouc/mite_extras)
- Run `hatch -v env create`

## Quick Start

### Run with `hatch`:
- `hatch run mite_extras`

## For devs

- Install developer dependencies with `hatch -v env create dev`
- Initialize `pre-commit` with `hatch run dev:pre-commit install`
- Run tests with `hatch run dev:pytest`
- Run CLI with `hatch run dev:mite_extas`
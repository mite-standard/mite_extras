mite_extras
==========

Validation and auxiliary functionality for MITE entries.

Attention: this program has only been tested with Ubuntu Linux.

## Installation

### With `hatch` from GitHub

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Download or clone the [repository](https://github.com/mmzdouc/mite_extras)
- Run `hatch -v env create`

## Quick Start: Example

### Run with `hatch`:

To validate MITE entries or update them to a new schema version

- `hatch run mite_extras -i input/ -o output/ -fout json`

## For devs

- Install developer dependencies with `hatch -v env create dev`
- Initialize `pre-commit` with `hatch run dev:pre-commit install`
- Run tests with `hatch run dev:pytest`
- Run CLI with `hatch run dev:mite_extras` and the appropriate options
- If necessary, remove the environment again with `hatch env remove dev`
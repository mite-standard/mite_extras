mite_extras
==========

Validation and conversion functionality for reaction SMARTS collected in the scope of MITE.

Attention: this program has only been tested with Ubuntu Linux.

## Installation

### With `hatch` from GitHub

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Download or clone the [repository](https://github.com/mmzdouc/mite_extras)
- Run `hatch -v env create`

## Quick Start: Example

### Run with `hatch`:

#### Input format: `raw` (MIBiG Submission Portal)

First, pre-format the raw input files:

- e.g. `hatch run python input_file_splitter.py example_input/new608.json`

Move the generated file(s) into a different directory (e.g. 'input')

- `hatch run mite_extras -i input/ -o output/ -fin raw -fout json`

This will convert all files into the MITE schema.

#### Input format: `mite` (MITE-formatted json)

To validate MITE entries or update them to a new schema version

- `hatch run mite_extras -i input/ -o output/ -fin mite -fout json`

## For devs

- Install developer dependencies with `hatch -v env create dev`
- Initialize `pre-commit` with `hatch run dev:pre-commit install`
- Run tests with `hatch run dev:pytest`
- Run CLI with `hatch run dev:mite_extras` and the appropriate options
- If necessary, remove the environment again with `hatch env remove dev`
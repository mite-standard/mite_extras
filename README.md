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

First, pre-format the raw input files (assuming that they are coming from the MIBiG Submission Portal):

- e.g. `hatch run python input_file_splitter.py example_input/new_608.json`

Move the generated file(s) into a different directory (e.g. 'input')

- `hatch run mite_extras -i input/ -o output/ -m raw`

This will convert all files into the MITE schema.

## For devs

- Install developer dependencies with `hatch -v env create dev`
- Initialize `pre-commit` with `hatch run dev:pre-commit install`
- Run tests with `hatch run dev:pytest`
- Run CLI with `hatch run dev:mite_extas`
mite_extras
==========

[![DOI](https://zenodo.org/badge/804997522.svg)](https://doi.org/10.5281/zenodo.13941745) 
[![PyPI version](https://badge.fury.io/py/mite-extras.svg)](https://badge.fury.io/py/mite-extras)


Contents
-----------------
- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Start](#quick-start)
- [Attribution](#attribution)
- [For Developers](#for-developers)

## Overview

**MITE** (Minimum Information about a Tailoring Enzyme) is a community-driven database for the characterization of tailoring enzymes.
These enzymes play crucial roles in the biosynthesis of secondary or specialized metabolites, naturally occurring molecules with strong biological activities, such as antibiotic properties.

This repository contains functionality to validate MITE entries in JSON format.

For more information, visit the [MITE Data Standard Organization page](https://github.com/mite-standard) or read our [publication]( https://doi.org/10.1093/nar/gkaf969).

## Documentation

This repository provides validation functionality for verifying files following the MITE `JSON Schema` format, such as the files in the  **MITE** data repository [mite_data](https://github.com/mite-standard/mite_data).

The validation functionality includes:
- Calling the `mite_schema` functionality
- Validation of reaction SMARTS (does the reaction SMARTS lead to the expected product when applied to a specified substrate)
- Sanitation of SMILES and reaction SMARTS
- Fetching of NCBI GenPept/UniProtKB IDs

`mite_extras` can be used as a CLI to automatically update MITE entries.
Furthermore, this repository can also be used as a library. For examples, see [mite_data](https://github.com/mite-standard/mite_data) or [mite_web](https://github.com/mite-standard/mite_web).

For errors, feature requests, and suggestions, please open an [Issue](https://github.com/mite-standard/mite_extras/issues) or start a discussion in the [MITE Discussion Forum](https://github.com/orgs/mite-standard/discussions/).

## System Requirements

### OS Requirements

Local installation was tested on:

- Ubuntu Linux 20.04 and 22.04 (command line)

#### Python dependencies

Dependencies including exact versions are specified in the [pyproject.toml](./pyproject.toml) file.

## Installation Guide

### With `pip`

- Install with `pip install mite_extras`

## Quick Start

To validate MITE entries or update them to a new schema versiom, run:

- `mite_extras -i <input/> -o <output/>`

## Attribution

### License

`mite_extras` is an open source tool licensed under the MIT license (see [LICENSE](LICENSE)).

### Publications

See [CITATION.cff](CITATION.cff) or [MITE online](https://mite.bioinformatics.nl/) for information on citing MITE.

### Acknowledgements

This work was supported by the Netherlands Organization for Scientific Research (NWO) KIC grant KICH1.LWV04.21.013.

## For Developers

*Nota bene: for details on how to contribute to the MITE project, please refer to [CONTRIBUTING](CONTRIBUTING.md).*

### Package Installation

*Please note that the development installation is only tested and supported on (Ubuntu) Linux.*

### With `uv` from GitHub

*Note: assumes that `uv` is installed locally - see the methods described [here](https://docs.astral.sh/uv/getting-started/installation/).* 

```commandline
git clone https://github.com/mite-standard/mite_extras
uv sync --extra dev
uv run pre-commit install
```

All tests should be passing
```commandline
uv run pytest
```

### CI/CD and Deployment

CI/CD via GitHub Actions runs on every PR and push to the `main` branch.

A new release created on the [mite_extras](https://github.com/mite-standard/mite_extras) GitHub page will automatically relay changes to [PyPI](https://pypi.org/project/mite-extras/) and [Zenodo](https://doi.org/10.5281/zenodo.13941744).

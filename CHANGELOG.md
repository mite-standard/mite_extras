# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.6.3] 30-09-2025

### Changed

- Dependency update
- Switched installer from `hatch` to `uv`
- `CITATION.cff` update to cover publication in *NAR*

## [1.6.2] 29-08-2025

### Changed

- `to_json()` and `to_html()` return changelog even if it is an empty list (compatibility with mite_web data submission)

## [1.6.1] 29-08-2025

### Changed

- Unpinned version of dependency of `mite_schema` to facilitate maintenance (should always pull the newest release)

## [1.6.0] 28-08-2025

### Changed

- In html_json generation, changed substrate from tuple to list of tuples to visualize substrates separately
- Changed Wikidata QID checking from raise errror to log warning (frequent 403 errors)

## [1.5.1] 26-07-2025

### Fixed

- Removed SMARTS atom index check due to backward incompatibility

## [1.5.0] 17-06-2025

### Added

- Compatibility with `mite_schema` version `1.8.0`
- Higher tolerance to malformed JSON files with `json-repair`
- Validation check for desired/undesired SMARTS patterns

### Changed

- Code cleanup

## [1.4.1] 13-02-2025

### Changed

- Bumped `mite_schema` version to `1.7.1`

## [1.4.0] 13-02-2025

### Changed

- Bumped `mite_schema` version to `1.7.0`
- Reworked `ValidationManager`
- Implemented option for intramolecular reactions (recognized from reaction SMARTS)
- Pinned RDkit version to `2024.3.6` 

## [1.3.1] 02-12-2024

### Added

- Fixed formatting of Ketcher-flavored charges in reaction SMARTS

### Bugfix

- Fixed UniProt fetching

## [1.3.0] 30-11-2024

### Changed

- Updated version of `mite_schema` to `1.6.0`
- Moved generation of databaseId hyperlinks from html_json() to `mite_web`

### Bugfix

- Fixed UniProt SPARQL fetching
- Fixed MIBiG URL formatting

## [1.2.1] 09-11-2024

### Changed

- Implemented additional checks to prevent writing of empty databaseID objects

## [1.2.0] 09-11-2024

### Changed

- Updated to mite_schema version 1.5.1

### Removed

- Removed backward compatibility to mite schema version <1.5.1

## [1.1.0] 21-10-2024

### Added

- Cross-link to Expansy Enzymes

### Changed

- Updated README to clarify scope of package

### Removed

- Removed id-mappings file: not needed anymore
- Removed html-generation option from CLI
- Removed html template generation
- Removed dependencies related to generation of static html pages

## [1.0.1] 16-10-2024

### Added

- Fixes to metadata of project

## [1.0.0] 16-10-2024

### Added

- Initial version of `mite_extras`
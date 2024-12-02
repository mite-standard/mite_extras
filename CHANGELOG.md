# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
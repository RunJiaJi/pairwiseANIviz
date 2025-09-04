# Changelog for pairwiseANIviz

## [1.2] - 2025-09-04

### Added
- Automatic header detection for FastANI output files
- Position-based column renaming for maximum compatibility
- **Universal whitespace delimiter support** using `\s+` regex pattern

### Changed
- Improved input file parsing to handle headers like:
  - `query reference ani mapped_fragments total_fragments`
  - `Query Reference ANI mapped_fragments total_fragments`
  - Any other format as long as first 3 columns are in correct order
- Enhanced error handling for malformed input files

### Fixed
- Function name typo in comments (readIput → readInput)
- Better handling of edge cases in header detection

## [1.1]

### Added
- Initial release of pairwiseANIviz
- Hierarchical clustering visualization of ANI results
- Support for various matplotlib colormaps
- Taxonomic classification integration with GTDB-Tk
- Multiple output formats (JPG, PNG, TIFF, SVG, PDF, EPS)
- Customizable clustering methods and distance metrics
- ANI value annotation and threshold highlighting

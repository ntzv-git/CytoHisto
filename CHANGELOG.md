# Changelog

All notable changes to CytoHisto. Versions follow [semantic versioning](https://semver.org).

## [1.1.1] - 2026-10-06
### Fixed
- Linux: figures were always saved as PNG, whatever the format selected in the save
  dialog, when the file name was typed without an extension. The extension now follows the
  selected format (PNG, TIFF, PDF, SVG); an extension typed in the name still takes priority.

## [1.1.0] - 2026-10-06
### Documentation
- README rewritten in Markdown (`README.md`) with screenshots; ready-to-use programs in `bin/`.

### Changed
- **Input channels are now set per file** (graphical application and R version).
  Each file is converted with its own input range, read from its `$PnR` keyword, so files
  recorded with different resolutions (e.g. 1024 and 65536 channels) are overlaid on the
  same output scale. Previously the range of the first loaded file was applied to all files.
- The file list shows the input range of each file.
- R version: `-i/--in-channels` is now a file option (placed after `-f`, or before the first
  `-f` to apply to every file).

### Added
- Version number: shown in the window title, `--version` option in both versions.

## [1.0.0] - 2026-10-05
First release.
- Graphical application (Windows, Linux, macOS) and R command-line version.
- FCS 2.0/3.0/3.1 reader, including parameters stored with different bit widths.
- Overlay of several files; per-curve label, colour, line width/type, smoothing, gate
  (%, mean, CV), parameter.
- Output channels, X zoom, automatic or manual axis steps, labels above peaks or legend.
- Image size in pixels, five independent text sizes scaling with the image width.
- Drag and drop, Ctrl+V paste, mouse coordinates, statistics (CSV), settings (JSON).
- PNG, TIFF, PDF and SVG export.

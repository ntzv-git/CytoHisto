# Changelog

All notable changes to CytoHisto. Versions follow [semantic versioning](https://semver.org).

## [1.3.2] - 2026-10-06
### Changed
- **Settings (JSON) store every imported file**, selected or not, with all its settings, plus
  the selection, which is restored when the settings are opened. (In 1.3.1 only the selected
  files were stored.) Figure and statistics still use the selected files only.
- Each file is stored with its absolute path and its path relative to the settings file.

### Added
- **Moved or deleted files no longer prevent opening settings.** Each file is looked for at its
  saved location, then relative to the settings file (a whole project folder can be moved),
  then next to it. Missing files can be searched for by name in a folder of your choice
  (sub-folders included); files still missing or unreadable are skipped and listed, and
  everything else is loaded.
- Clear error messages for files that are not FCS files or are truncated.

## [1.3.1] - 2026-10-06
### Changed
- The file list is taller (8 rows) and has a scroll bar; the mouse wheel scrolls the list
  when the pointer is over it (the settings panel otherwise).
- Every save uses the selected files only: figure, statistics (CSV) and settings (JSON).

## [1.3.0] - 2026-10-06
### Added
- **Only the selected files are plotted.** The file list accepts multiple selection
  (Ctrl+click, Shift+click, Ctrl+A); the plot, the statistics and the saved figure follow the
  selection. Newly imported files join the selection. The Y axis is recalculated on the
  highest visible peak at every change of selection, unless a Y zoom is set.
- "Remove" removes every selected file; "Up"/"Down" move the active file (the one shown in
  "Selected curve", i.e. the last one clicked).

## [1.2.0] - 2026-10-06
### Changed
- **Histogram resolution is now fixed at 1024 classes over the full scale**, whatever the
  output channels, which only set the X axis unit. Changing the output channels no longer
  changes the Y axis, and high resolutions (e.g. 65536 channels) no longer turn the curves
  into unreadable, overlapping noise. Smoothing now merges classes (1 = 1024 classes).
- Output channels are chosen from a list; when they change, the X zoom, X steps and gates
  are rescaled so that they keep the same place on the curves.

### Added
- **Y zoom** (minimum and maximum of the Y axis); `-Y/--ylim` in the R version.

### Fixed
- The X zoom can no longer go beyond the channel range (0 to output channels): it is limited
  and the field is corrected.

## [1.1.2] - 2026-10-06
### Fixed
- Standalone programs (Linux and Windows builds): saving as PDF or SVG failed with
  "No module named 'matplotlib.backends.backend_pdf'" (or `backend_svg`). These writers are
  now imported explicitly so that PyInstaller bundles them. PNG and TIFF were not affected,
  nor was `python cytohisto.py`.

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

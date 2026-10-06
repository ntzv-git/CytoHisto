CytoHisto - publication-ready flow cytometry histograms
=======================================================
Version 1.1.0 (see CHANGELOG.md). "python cytohisto.py --version" / "Rscript cytohisto.R --version".

Works with any FCS 2.0 / 3.0 / 3.1 list-mode file, including files whose parameters are
stored with different bit widths (which flowCore cannot read).

Files
  cytohisto.py            graphical application (Windows, Linux, macOS)
  build_windows_exe.bat   builds CytoHisto.exe on Windows (once)
  build_linux.sh          builds a standalone CytoHisto program on Linux (once)
  cytohisto.R             command-line version for R (needs the ggplot2 package)

Build the Windows executable (once, ~5 min)
  1. Install Python for Windows: https://www.python.org/downloads/
     -> during installation, tick "Add python.exe to PATH".
  2. Double-click build_windows_exe.bat.
  3. The program is dist\CytoHisto.exe. It runs on its own and can be copied to other
     Windows PCs without Python.
  On first launch, Windows SmartScreen may show "Windows protected your PC" (unsigned program):
  "More info" -> "Run anyway".

Build the Linux program (once)
  sudo apt install python3-tk python3-venv
  ./build_linux.sh            -> dist/CytoHisto (runs on its own)
  It uses the system Python on purpose: the Tk library shipped with conda has no
  anti-aliased fonts, which makes the interface look rough.

Run without building (any system)
  pip install numpy matplotlib tkinterdnd2
  python cytohisto.py
  (on Linux, also install the python3-tk package; prefer the system Python to conda,
  or in conda run: conda install -c conda-forge "tk=*=xft_*")

Graphical application
  - Plot settings (all curves): axis titles, output channels (X axis resolution), X zoom, axis steps (empty = automatic; e.g.
    "200,100,20" for X = labels, grid, ticks and "100,20" for Y = labels, ticks; steps
    too small for the range are ignored), names above the peaks or legend.
  - Image: size in pixels (default 2000 x 1000) and five text sizes (labels, X title,
    Y title, X values, Y values); text, lines and ticks scale with the image width.
  - FCS files: drag and drop .fcs files anywhere in the window, click "Add...", or copy
    the files in the file manager and press Ctrl+V in the window.
    Drops are reported in the terminal ("[CytoHisto] drag entered" / "drop received").
    They are overlaid; the Y axis fits the highest visible peak of all files.
  - The file list shows each curve's colour, file name and label.
  - Selected curve: click a file in the list to set its label (name on the plot and in the
    statistics), colour (#D55E00, #d50,
    rgb(213,94,0), crimson... or "Pick..."), line width, line type, smoothing, gate
    (highlighted region + %, mean, CV), parameter and input channels.
  - Input channels are set per file: each file is read with its own range ($PnR, shown in
    the list), so files from instruments with different resolutions (e.g. 1024 and 65536)
    are converted to the same output channels. The value can be changed if a file declares
    a wrong range.
  - The preview updates automatically, with the exact proportions of the saved image.
  - Move the mouse over the plot: a cursor line follows it and the line above the plot shows
    the cursor position (x, y).
  - Save figure (PNG, TIFF, PDF, SVG), statistics (CSV) and all settings (JSON).

R command-line version
  Rscript cytohisto.R -h        lists all options. Example:
  Rscript cytohisto.R -L -X 50,250 -s 2 -o comparison.png -P 3000x1500 \
      -f 225.fcs -n "CIRAD225" -c "#0072B2" -g 195,235 \
      -f 313.fcs -n "CIRAD313" -c crimson -w 0.8 -g 85,105

#!/usr/bin/env bash
# Builds a standalone CytoHisto program for Linux (dist/CytoHisto). Run it once, from this folder.
# Uses the system Python on purpose: conda's Tk has no anti-aliased fonts (ugly interface).
# Requirement (Debian/Ubuntu): sudo apt install python3-tk python3-venv
set -e
cd "$(dirname "$0")"
PY="${PYTHON:-/usr/bin/python3}"
"$PY" -c "import tkinter" 2>/dev/null || {
  echo "Tkinter is missing for $PY. Install it with: sudo apt install python3-tk python3-venv"; exit 1; }

echo "=== Creating an isolated environment (.venv-linux) with $PY..."
"$PY" -m venv .venv-linux
. .venv-linux/bin/activate

echo "=== Installing numpy, matplotlib, tkinterdnd2 and pyinstaller..."
python -m pip install --upgrade pip >/dev/null
python -m pip install numpy matplotlib pillow tkinterdnd2 pyinstaller

echo "=== Building the program..."
pyinstaller --noconfirm --clean --onefile --windowed --name CytoHisto --collect-all tkinterdnd2 cytohisto.py

echo
echo "Done: dist/CytoHisto  (runs on its own; double-click it or run ./dist/CytoHisto)"

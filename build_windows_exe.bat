@echo off
REM Builds CytoHisto.exe (Windows). Run it once, from this folder.
REM Requirement: Python 3.10+ for Windows (https://www.python.org, tick "Add python.exe to PATH").

cd /d "%~dp0"
where python >nul 2>nul || (echo Python was not found. Install it from https://www.python.org & pause & exit /b 1)

echo === Creating an isolated environment...
python -m venv .venv || (pause & exit /b 1)
call .venv\Scripts\activate.bat

echo === Installing numpy, matplotlib, tkinterdnd2 and pyinstaller...
python -m pip install --upgrade pip >nul
python -m pip install numpy matplotlib pillow tkinterdnd2 pyinstaller || (pause & exit /b 1)

echo === Building the executable...
pyinstaller --noconfirm --clean --onefile --windowed --name CytoHisto --collect-all tkinterdnd2 cytohisto.py || (pause & exit /b 1)

echo.
echo Done: dist\CytoHisto.exe
echo This .exe runs on its own: copy it anywhere (Python is no longer needed).
explorer dist
pause

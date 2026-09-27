@echo off
setlocal enabledelayedexpansion

REM Set directory context to the location of this batch script
set SCRIPT_DIR=%~dp0

echo ==============================================================================
echo Building Vivado Project using range_detector.tcl...
echo ==============================================================================

call settings64.bat
@echo Vivado TCL script is running...
@echo Check log files for progress.
call vivado -nolog -nojournal -mode batch -source "%SCRIPT_DIR%range_detector.tcl" -tclargs %1 > "%SCRIPT_DIR%range_detector.log"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Vivado project generation failed. Check range_detector.log for details.
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo [SUCCESS] Project created successfully!
echo.

set /p OPEN_GUI="Do you want to open the project in Vivado GUI now? (Y/N): "
if /i "%OPEN_GUI%"=="Y" (
    echo Opening Vivado GUI...
    start "" vivado "%SCRIPT_DIR%..\project\range_detector_proj.xpr"
)

endlocal
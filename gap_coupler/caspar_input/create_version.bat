@echo off
setlocal enabledelayedexpansion

:: Navigate three directories up
set "VERSION_FILE=%~dp0..\..\..\version.json"

:: Read the version from version.json
for /f "usebackq tokens=2 delims=:," %%A in ("%VERSION_FILE%") do (
    set VERSION=%%A
    set VERSION=!VERSION:~2,-1!
    
    :: Create the version_define.h file in the same directory as the script
    echo #define PROJECT_VERSION "!VERSION!" > "%~dp0version_define.h"
    echo Created version_define.h with PROJECT_VERSION=!VERSION!
)

endlocal

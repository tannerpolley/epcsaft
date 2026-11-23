@echo on
call "%RECIPE_DIR%\build.bat"
if errorlevel 1 exit /b 1
"%PYTHON%" "%RECIPE_DIR%\qualification\retain.py" %JOBS%
if errorlevel 1 exit /b 1

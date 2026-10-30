@echo on
rem Qualification build: the shipped build plus the qualification programs and development records.
%PYTHON% packaging/conda-forge/qualification/prepare.py
if errorlevel 1 exit /b 1
call "%RECIPE_DIR%\build.bat"
if errorlevel 1 exit /b 1
setlocal
call run_autotools_clang_conda_build.bat packaging/conda-forge/qualification/hs071.sh
if errorlevel 1 exit /b 1
endlocal
%PYTHON% packaging/conda-forge/qualification/retain.py
if errorlevel 1 exit /b 1

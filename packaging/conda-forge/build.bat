@echo on
set EPCSAFT_BUILD_JOBS=6
setlocal
call run_autotools_clang_conda_build.bat packaging/conda-forge/qualification/build-ipopt-only.sh
if errorlevel 1 exit /b 1
endlocal
set CC=clang-cl
set CXX=clang-cl
set CFLAGS=/MD /fp:precise /clang:-ffp-contract=off
set CXXFLAGS=/MD /fp:precise /clang:-ffp-contract=off
set PKG_CONFIG_PATH=%LIBRARY_PREFIX%\lib\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig
set PATH=%PREFIX%\bin;%LIBRARY_PREFIX%\bin;%PATH%
%PYTHON% packaging/conda-forge/qualification/prepare.py
if errorlevel 1 exit /b 1
%PYTHON% -m pip install ./engine -vv --no-deps --no-build-isolation --config-settings=build-dir=%SRC_DIR%\engine\build-win --config-settings=setup-args=--wrap-mode=nodownload --config-settings=compile-args=-j6
if errorlevel 1 exit /b 1
copy %LIBRARY_PREFIX%\bin\epcsaft_ipopt_3119.dll %SP_DIR%\epcsaft\
%PYTHON% packaging/conda-forge/qualification/retain.py
if errorlevel 1 exit /b 1

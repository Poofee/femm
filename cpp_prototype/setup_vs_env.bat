@echo off
echo Setting up Visual Studio 2022 Developer Command Prompt...

REM 设置Visual Studio环境变量
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"

echo Environment setup complete.
echo Available compilers:
where cl

echo.
echo Now you can compile the test programs:
echo cl /EHsc minimal_test.cpp
echo cl /EHsc verify_core.cpp

pause
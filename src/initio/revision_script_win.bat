@echo off
rem Windows batch file for retrieving svn version information
SETLOCAL
for /f "delims=*" %%V in ('svnversion .\') do set REV=%%V
echo #define ADHREV "%REV%" > %1
for /F "usebackq tokens=1,2 delims==" %%i in (`wmic os get LocalDateTime /VALUE 2^>NUL`) do if '.%%i.'=='.LocalDateTime.' set ldt=%%j
echo #define ADHREVDATE "%ldt:~0,4%.%ldt:~4,2%.%ldt:~6,2%" >> %1
echo #define ADHREVTIME "%ldt:~8,2%:%ldt:~10,2%:%ldt:~12,2%" >> %1
exit 0
ENDLOCAL
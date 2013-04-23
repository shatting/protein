@echo off

echo.---------------------------------------------------------------
echo.------------------------  RUNNING LGA -------------------------
echo.---------------------------------------------------------------

".\tools\winscp\winscp.com" /script:"scripts\runlga.txt" /privatekey="c:\login.ppk"
@echo off

".\tools\winscp\winscp.com" /script:"#{lscriptdir}#\sync.txt" /privatekey="#{keyfile}#"
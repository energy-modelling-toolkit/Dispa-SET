@echo off

IF EXIST gams_path.txt SET /p GAMSPATH=< gams_path.txt
if not exist gams_path.txt SET GAMSPATH=C:\GAMS\win32\24.3
:startloop
set msg=Provided gams path not found: 
set msg=%msg%%GAMSPATH%
IF NOT EXIST %GAMSPATH% (
	echo %msg%
	echo %GAMSPATH%>gams_path.txt
	SET /p GAMSPATH="Enter Gams path (e.g. C:\GAMS\win32\24.3):"
	pause
	goto startloop
)

IF EXIST %GAMSPATH% echo %GAMSPATH%>gams_path.txt
IF EXIST %GAMSPATH% SET fullpath=%gamspath%\gams.exe
if exist Inputs.gdx del Inputs.gdx

%fullpath% make_gdx.gms

if not exist Inputs.gdx (
	del *.gdx
	IF EXIST make_gdx.lst del make_gdx.lst
	if exist make_gdx.log del make_gdx.log
	if exist make_gdx.lxi del make_gdx.lxi
	echo.
	echo Something has gone wrong. Run make_gdx.gms in Gams to see the error.
	pause
)

if exist Inputs.gdx (
	rename Inputs.gdx Inputs
	del *.gdx
	rename Inputs Inputs.gdx
	IF EXIST make_gdx.lst del make_gdx.lst
	if exist make_gdx.log del make_gdx.log
	if exist make_gdx.lxi del make_gdx.lxi
	echo.
	echo DispaSET inputs successfully written to Inputs.gdx
	pause
)

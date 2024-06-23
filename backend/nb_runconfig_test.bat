@ECHO OFF

REM This script records some details about your NumBAT environment
REM to help identify the cause of any crashes or other difficulties

SET OUT=.\nb_build_config.txt

IF EXIST "backend\" (
  ECHO This script must be executed from the NumBAT root directory.
  ECHO.
  EXIT /B 1
)

ECHO.> "%OUT%"

ECHO NumBAT config test file >> "%OUT%"
ECHO.-------------- >> "%OUT%"
ECHO. >> "%OUT%"
ECHO Date:        >> "%OUT%" DATE /T >> "%OUT%" 
ECHO Env:        >> "%OUT%" WMIC OS Get /Format:List >> "%OUT%"

ECHO python path:  >> "%OUT%" WHERE python >> "%OUT%"
ECHO python3 path:  >> "%OUT%" WHERE python >> "%OUT%"
ECHO python3 ver.: >> "%OUT%" python -c "import sys; print(sys.version)" >> "%OUT%"
ECHO numpy ver.:   >> "%OUT%" python -c "import numpy; print(numpy.__version__)" >> "%OUT%"
ECHO scipy ver.:   >> "%OUT%" python -c "import scipy; print(scipy.__version__)" >> "%OUT%"
ECHO matplotlib ver.: >> "%OUT%" python -c "import matplotlib; print(matplotlib.__version__)" >> "%OUT%"

ECHO. >> "%OUT%"

ECHO. >> "%OUT%"
ECHO NumBAT build settings: >> "%OUT%"
ECHO. >> "%OUT%"
TYPE ./backend/fortran/nb_build_config.txt >> "%OUT%"

ECHO Done! Information written to "%OUT%".


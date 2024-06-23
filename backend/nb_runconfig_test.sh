#!/bin/sh

# This program records some details about your NumBAT environment
# to help identify the cause of any crashes or other difficulties
#

OUT="./nb_config.txt"

if [ ! -e "backend" ]; then
    echo
    echo This script must be executed from the NumBAT root directory.
    echo
    exit
fi

echo > ${OUT}

echo NumBAT config test file  >> ${OUT}
echo ------------------------  >> ${OUT}
echo >> ${OUT}
echo "Date:               " `date` >> ${OUT}
echo "Env:                " ` uname -a` >> ${OUT}
echo "lib search path:    " `set|grep LD_LIBRARY_PATH` >> ${OUT}

echo "python path:        " `which python` >> ${OUT}
echo "python3 path:       " `which python` >> ${OUT}
echo "python3 ver.:       " `python -c "import sys;print(sys.version)" ` >> ${OUT}
echo "numpy ver.:         " `python -c "import numpy;print(numpy.__version__)" ` >> ${OUT}
echo "scipy ver.:         " `python -c "import scipy;print(scipy.__version__)" ` >> ${OUT}
echo "matplotlib ver.:    " `python -c "import matplotlib;print(matplotlib.__version__)" ` >> ${OUT}

echo "core dumps enabled: " `ulimit -c` >>${OUT}
echo "core dump location: " `cat /proc/sys/kernel/core_pattern ` >>${OUT}

echo >> ${OUT}
echo "nb_fortran.so linking" >> ${OUT}
ldd backend/fortran/nb_fortran.so >> ${OUT}

echo >> ${OUT}
echo "NumBAT build settings: " >> ${OUT}
echo >> ${OUT}
cat ./backend/fortran/nb_build_config.txt >> ${OUT}



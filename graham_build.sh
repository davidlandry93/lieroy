#!/usr/bin/env python3

cmake .. -DCMAKE_PREFIX_PATH=~/local -DCMAKE_INSTALL_PREFIX=~/local -DBOOST_ROOT=~/local -DBUILD_PYTHON_BINDINGS=1 -DPYTHON_INCLUDE_OVERRIDE=/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/include/python3.5m


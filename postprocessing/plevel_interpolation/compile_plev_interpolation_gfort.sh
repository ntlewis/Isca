#!/usr/bin/env bash
#Script that compiles the plev.x executable

cd ./exec

source $GFDL_BASE/src/extra/env/$GFDL_ENV_PP

../bin/mkmf -p plev.x -t ../bin/mkmf.template.gfort -c "-Duse_netCDF" -a ../src ../src/path_names ../src/shared/mpp/include ../src/shared/include

make -f Makefile

#!/bin/sh
###################################################################
##                                                               ##
##              c  o  n  f  i  g  u  r  e  .  s  h               ##
##                                                               ##
## configure script to link to appropriate Makefile              ##
##                                                               ##
## written by Werner von Bloh                                    ##
## Potsdam Institute for Climate Impact Research                 ##
## P.O. Box 60 12 03                                             ##
## 14412 Potsdam, Germany                                        ##
##                                                               ##
## Last change: 22.10.2004                                       ##
##                                                               ##
###################################################################

echo Configuring lpj...
        
osname=$(uname)
if [ "$osname" = "Linux" ] 
then
  echo Operating system is Linux
  ln -sf Makefile.gcc Makefile.inc
  exit 0
fi
if [ "$osname" = "AIX" ] 
then
  echo  Operating system is AIX
  ln -sf Makefile.aix Makefile.inc
  exit 0
fi
echo >&2 Warning: unsupported operating system, Makefile.$osname created
cp Makefile.gcc Makefile.$osname
ln -sf Makefile.$osname Makefile.inc
exit 1

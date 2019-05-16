########################################################################
##                                                                    ##
##               M  a  k  e  f  i  l  e                               ##
##                                                                    ##
##   Makefile for C implementation of LPJ                             ##
##                                                                    ##
##   written by Werner von Bloh                                       ##
##   Potsdam Institute for Climate Impact Research                    ##
##   P.O. Box 60 12 03                                                ##
##   14412 Potsdam/Germany                                            ##
##                                                                    ##
##   Last change: 19.11.2004                                          ##
##                                                                    ##
########################################################################

include Makefile.inc

INC     = include

HDRS    = $(INC)/buffer.h $(INC)/coord.h $(INC)/lpj.h $(INC)/pftlist.h\
          $(INC)/soil.h $(INC)/climate.h $(INC)/date.h $(INC)/pft.h\
          $(INC)/pftpar.h $(INC)/types.h $(INC)/tree.h $(INC)/header.h\
          $(INC)/crop.h $(INC)/grass.h $(INC)/errmsg.h $(INC)/numeric.h\
          $(INC)/conf.h $(INC)/swap.h $(INC)/soilpar.h $(INC)/stand.h\
          $(INC)/list.h $(INC)/cell.h  $(INC)/units.h $(INC)/output.h\
          $(INC)/config.h $(INC)/popen.h

CONF	= lpj.conf output.conf cru_1901-2003.conf

DATA    = par/soil.par par/pft.par

CRU	= input/soil.bin input/co2_2003.dat input/cld2003.clm\
          input/tmp2003.clm input/pre2003.clm\
          input/wet2003.clm input/grid.bin

SCRIPTS	= distribute.sh configure.sh newpft.sh

FILES	= $(HDRS) Makefile Makefile.aix Makefile.gcc Makefile.win32\
          $(CONF) $(DATA) $(SCRIPTS) INSTALL AUTHORS VERSION configure.bat

bin: 
	$(MKDIR) lib
	(cd src && $(MAKE))
cat:
	(cd src && $(MAKE) ../catlpj)
clm:
	(cd src && $(MAKE) ../cru2clm)
print:
	(cd src && $(MAKE) ../printlpj)
	

all: bin cat clm print


test: bin
	$(MKDIR) output
	lpj

clean:
	(cd src  && $(MAKE) clean)

archive: 
	tar -cf LPJ-Cv3.tar $(FILES) src/base/*.c src/base/Makefile\
        src/tree/*.c src/soil/*.c src/soil/Makefile src/tree/Makefile\
        src/grass/*.c src/grass/Makefile src/numeric/*.c src/numeric/Makefile\
        src/Makefile src/*.c  src/lpj/*.c src/lpj/Makefile\
        src/climate/Makefile src/climate/*.c
	gzip -f LPJ-Cv3.tar

zipfile: 
	zip -l LPJ-Cv3.zip $(FILES) src/base/*.c src/base/Makefile\
        src/tree/*.c src/soil/*.c src/soil/Makefile src/tree/Makefile\
        src/grass/*.c src/grass/Makefile src/numeric/*.c src/numeric/Makefile\
        src/Makefile src/*.c  src/lpj/*.c src/lpj/Makefile\
        src/climate/Makefile src/climate/*.c

cru.zip: $(CRU)
	zip cru.zip $(CRU)

cru.tar: $(CRU)
	tar -cf cru.tar $(CRU)

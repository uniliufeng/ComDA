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

include ../Makefile.inc

OBJS	= cru2clm.$O printlpj.$O catlpj.$O lpj.$O

SRC	= cru2clm.c printlpj.c catlpj.c

INC     = ../include

HDRS    = $(INC)/buffer.h $(INC)/coord.h $(INC)/lpj.h $(INC)/pftlist.h\
          $(INC)/soil.h $(INC)/climate.h $(INC)/date.h $(INC)/pft.h\
          $(INC)/pftpar.h $(INC)/types.h $(INC)/tree.h $(INC)/header.h\
          $(INC)/crop.h $(INC)/grass.h $(INC)/errmsg.h $(INC)/numeric.h\
          $(INC)/conf.h $(INC)/swap.h $(INC)/soilpar.h $(INC)/stand.h\
          $(INC)/list.h $(INC)/cell.h  $(INC)/units.h $(INC)/output.h\
          $(INC)/config.h $(INC)/popen.h

LIB	= ../lib

LPJLIBS	= $(LIB)/liblpj.$A  $(LIB)/libbase.$A  $(LIB)/libsoil.$A\
          $(LIB)/libtree.$A  $(LIB)/libgrass.$A  $(LIB)/libnum.$A\
          $(LIB)/libclimate.$A

EXE	= ../lpj$E

bin: 
	(cd base && $(MAKE))
	(cd climate && $(MAKE))
	(cd numeric && $(MAKE))
	(cd lpj && $(MAKE))
	(cd soil && $(MAKE))
	(cd tree && $(MAKE))
	(cd grass && $(MAKE))
	$(MAKE) lpj
clean:
	(cd base && $(MAKE) clean)
	(cd climate && $(MAKE) clean)
	(cd lpj && $(MAKE) clean)
	(cd numeric && $(MAKE) clean)
	(cd soil && $(MAKE) clean)
	(cd tree && $(MAKE) clean)
	(cd grass && $(MAKE) clean)
	$(RM) $(RMFLAGS) $(OBJS) ../lpj$E ../catlpj$E ../printlpj$E

$(OBJS): $(HDRS)

.c.$O: 
	$(CC) $(CFLAGS) -I$(INC) -c $*.c

lpj: $(LPJLIBS) lpj.$O
	$(LINK) $(LNOPTS)$(EXE) lpj.$O $(LPJLIBS)  $(LIBS)  

../printlpj: printlpj.$O $(LPJLIBS)
	$(LINK) $(LNOPTS) ../printlpj$E printlpj.$O $(LIBS) $(LPJLIBS)

../cru2clm$E: cru2clm.$O $(INC)/types.h $(INC)/swap.h $(INC)/climate.h\
        $(INC)/date.h $(INC)/buffer.h $(INC)/coord.h $(INC)/config.h\
        $(LIB)/libnum.$A $(LIB)/libclimate.$A
	$(LINK) $(LNOPTS)../cru2clm$E cru2clm.$O $(LIB)/libclimate.$A  $(LIB)/libnum.$A 

../catlpj$E: catlpj.$O $(INC)/types.h $(INC)/swap.h $(INC)/header.h\
             $(LIB)/libnum.$(A)
	$(LINK) $(LNOPTS)../catlpj$E catlpj.$(O) $(LIB)/libnum.$(A) 

all: $(EXE) ../printlpj ../cru2clm ../catlpj

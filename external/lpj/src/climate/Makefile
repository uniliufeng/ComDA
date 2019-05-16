########################################################################
##                                                                    ##
##               M  a  k  e  f  i  l  e                               ##
##                                                                    ##
##   Makefile for C implementation of climate part of LPJ             ##
##                                                                    ##
##   written by Werner von Bloh                                       ##
##   Potsdam Institute for Climate Impact Research                    ##
##   P.O. Box 60 12 03                                                ##
##   14412 Potsdam/Germany                                            ##
##                                                                    ##
##   Last change: 26.05.2005                                          ##
##                                                                    ##
########################################################################

include ../../Makefile.inc

OBJ	= initclimate.$O getclimate.$O freeclimate.$O\
          prdaily.$O getclimate.$O coord.$O climateheader.$O getco2.$O

INC	= ../../include
LIB	= ../../lib/libclimate.$A

HDRS	= $(INC)/buffer.h $(INC)/coord.h $(INC)/lpj.h\
          $(INC)/climate.h $(INC)/date.h\
          $(INC)/pftpar.h $(INC)/types.h $(INC)/header.h\
          $(INC)/errmsg.h $(INC)/numeric.h\
          $(INC)/conf.h $(INC)/swap.h $(INC)/soilpar.h \
          $(INC)/list.h $(INC)/cell.h  $(INC)/units.h\
          $(INC)/config.h


$(LIB): $(OBJ)
	$(AR) $(ARFLAGS)$(LIB) $(OBJ)

$(OBJ): $(HDRS)

.c.$O: 
	$(CC) $(CFLAGS) -I$(INC) -c $*.c

clean: 
	$(RM) $(RMFLAGS) $(OBJ) $(LIB)

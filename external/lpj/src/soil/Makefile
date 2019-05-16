########################################################################
##                                                                    ##
##               M  a  k  e  f  i  l  e                               ##
##                                                                    ##
##   Makefile for C implementation of soil part of LPJ                ##
##                                                                    ##
##   written by Werner von Bloh                                       ##
##   Potsdam Institute for Climate Impact Research                    ##
##   P.O. Box 60 12 03                                                ##
##   14412 Potsdam/Germany                                            ##
##                                                                    ##
##   Last change: 05.11.2004                                          ##
##                                                                    ##
########################################################################

include ../../Makefile.inc

OBJ	= initsoil.$O soil.$O waterbalance.$O equilsoil.$O littersom.$O\
          getlag.$O temp_response.$O fscansoilpar.$O soiltemp.$O snow.$O\
          fire_sum.$O fire_prob.$O freadsoil.$O fwritesoil.$O fprintsoil.$O

INC     = ../../include
LIB	= ../../lib/libsoil.$A

HDRS    = $(INC)/buffer.h $(INC)/coord.h $(INC)/lpj.h $(INC)/pftlist.h\
          $(INC)/soil.h $(INC)/climate.h $(INC)/date.h $(INC)/pft.h\
          $(INC)/pftpar.h $(INC)/types.h\
          $(INC)/crop.h $(INC)/errmsg.h $(INC)/numeric.h\
          $(INC)/conf.h $(INC)/swap.h $(INC)/soilpar.h $(INC)/stand.h\
          $(INC)/list.h $(INC)/cell.h  $(INC)/units.h $(INC)/output.h\
          $(INC)/config.h $(INC)/popen.h

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS)$(LIB) $(OBJ)

$(OBJ): $(HDRS)

.c.$O: 
	$(CC) $(CFLAGS) -I$(INC) -c $*.c

clean: 
	$(RM) $(RMFLAGS) $(OBJ) $(LIB)  

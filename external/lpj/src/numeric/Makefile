########################################################################
##                                                                    ##
##               M  a  k  e  f  i  l  e                               ##
##                                                                    ##
##   Makefile for C implementation of numeric part of LPJ             ##
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

OBJ	= leftmostzero.$O bisect.$O linreg.$O fail.$O date.$O interpolate.$O\
          swap.$O list.$O buffer.$O header.$O fscanreal.$O fscanstring.$O\
          rand.$O popen.$O petpar.$O

INC	= ../../include
LIB	= ../../lib/libnum.$A

HDRS	= $(INC)/types.h $(INC)/numeric.h $(INC)/date.h $(INC)/swap.h\
          $(INC)/list.h $(INC)/buffer.h $(INC)/errmsg.h $(INC)/header.h\
          $(INC)/popen.h

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS)$(LIB) $(OBJ)

$(OBJ): $(HDRS)

.c.$O: 
	$(CC) $(CFLAGS) -I$(INC) -c $*.c

clean: 
	$(RM) $(RMFLAGS) $(OBJ) $(LIB)

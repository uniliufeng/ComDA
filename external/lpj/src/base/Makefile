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
##   Last change: 18.11.2004                                          ##
##                                                                    ##
########################################################################

include ../../Makefile.inc

OBJ	= fscanpftpar.$O firepft.$O climbuf.$O pftlist.$O\
          newpft.$O survive.$O establish.$O\
          fprintstand.$O temp_stress.$O\
          updategdd.$O fwritecell.$O fwriterestart.$O freadcell.$O\
          fscanconfig.$O nomortality.$O noinit.$O nofire.$O\
          freepft.$O photosynthesis.$O freegrid.$O killpft.$O\
          light.$O gp.$O gp_sum.$O water_stressed.$O interception.$O fpc_sum.$O\
          establishmentpft.$O freadstandlist.$O noadjust.$O\
          freadpft.$O fwritepft.$O fwritestand.$O\
          fprintpft.$O fprintstandlist.$O\
          fprintcell.$O writecoords.$O fwritestandlist.$O\
          freadstand.$O equilsom.$O initgdd.$O initpft.$O

INC	= ../../include
LIB	= ../../lib/libbase.$A

HDRS	= $(INC)/buffer.h $(INC)/coord.h $(INC)/lpj.h $(INC)/pftlist.h\
          $(INC)/soil.h $(INC)/climate.h $(INC)/date.h $(INC)/pft.h\
          $(INC)/pftpar.h $(INC)/types.h $(INC)/tree.h $(INC)/header.h\
          $(INC)/errmsg.h $(INC)/numeric.h\
          $(INC)/conf.h $(INC)/swap.h $(INC)/soilpar.h $(INC)/stand.h\
          $(INC)/list.h $(INC)/cell.h  $(INC)/units.h \
          $(INC)/config.h $(INC)/popen.h


$(LIB): $(OBJ)
	$(AR) $(ARFLAGS)$(LIB) $(OBJ)

$(OBJ): $(HDRS)

.c.$O: 
	$(CC) $(CFLAGS) -I$(INC) -c $*.c

clean: 
	$(RM) $(RMFLAGS) $(OBJ) $(LIB)

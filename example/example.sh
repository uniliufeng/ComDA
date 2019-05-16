#!/bin/sh

#lorenz
g++ -I../include lorenz.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/util/exception.cpp  -o test
g++ -I../include lorenz4d.cpp ../src/model/lorenz4d.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/util/exception.cpp  -o test -larmadillo

## enkf with lorenz
g++ -I../include enkf_lorenz.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/enkf.cpp  -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp -o test
# with armadillo lib
g++ -I../include enkf_lorenz_arma.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/ensemblekalmanfilter.cpp ../src/da/filter.cpp -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp -larmadillo -o test

# ukf with lorenz
g++ -I../include ukf_lorenz.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/unscentedkalmanfilter.cpp ../src/da/filter.cpp -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp -larmadillo -o test

# upf with lorenz
g++ -I../include upf_lorenz.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/unscentedparticlefilter.cpp ../src/da/particlefilter.cpp  -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp  ../src/da/filter.cpp -larmadillo -o test

# pf with lorenz
g++ -I../include pf_lorenz.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/particlefilter.cpp  -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp ../src/da/filter.cpp -larmadillo -o test
g++ -I../include pf_lorenz1.cpp ../src/model/lorenz.cpp ../src/model/model.cpp ../src/util/time.cpp ../src/da/particlefilter.cpp  -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/util/exception.cpp ../src/util/random.cpp ../src/util/rng.cpp ../src/da/filter_itpp.cpp -litpp -o test

# Geotop1.0
rm *.o
gcc -I../external/geotop/GEOtop/ -I ../external/geotop/Libraries/ASCII/ -I ../external/geotop/Libraries/FLUIDTURTLES/ -I ../external/geotop/Libraries/KeyPalette/ -I ../external/geotop/Libraries/GEOMORPHOLOGYLIB/ -I ../external/geotop/Libraries/MATH/   ../external/geotop/Libraries/ASCII/*.c ../external/geotop/Libraries/FLUIDTURTLES/*.c ../external/geotop/Libraries/KeyPalette/*.c ../external/geotop/Libraries/GEOMORPHOLOGYLIB/*.c ../external/geotop/Libraries/MATH/*.c ../external/geotop/GEOtop/*.c -c
rm geotop.o
g++ -I../external/geotop/GEOtop/ -I ../external/geotop/Libraries/ASCII/ -I ../external/geotop/Libraries/FLUIDTURTLES/ -I ../external/geotop/Libraries/KeyPalette/ -I ../external/geotop/Libraries/GEOMORPHOLOGYLIB/ -I ../external/geotop/Libraries/MATH/   *.o geotop.cpp -I../include ../src/model/geotop.cpp ../src/model/model.cpp ../src/util/time.cpp -o geotop/geotop -lm

# Vic-Nl
rm *.o
gcc -c -I ../external/vic/ ../external/vic/*.c
rm vicNl.o
g++ -I../include -I../external/vic vic.cpp ../src/model/vic.cpp ../src/model/model.cpp ../src/util/time.cpp *.o -o vic/vic

# Vic-Nl C++ version
# --no-warning
g++ -I../include -I../external/vic_c++ ../src/model/vic.cpp vic.cpp ../src/model/model.cpp ../src/util/time.cpp ../external/vic_c++/*.cpp -o vic/vic

# sib2
g++ -I../include sib2.cpp ../src/model/sib2model.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/model/model.cpp -o sib2/sib2

# colm
g++ -I../include colm.cpp ../src/model/commonlandmodel.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/util/landgrid.cpp ../src/model/model.cpp -o colmtest

# shaw
rm *.o
gfortran ../external/shaw/shaw23.for -c
g++ -I../include shaw.cpp ../src/model/shaw.cpp ../src/model/model.cpp ../src/util/time.cpp shaw23.o ../src/util/exception.cpp -o shaw/shaw -lgfortran

# lpj
# need make LPJ lib first
g++ -I../include -I../external ../src/model/lpj.cpp ../src/model/model.cpp ../src/util/exception.cpp ../src/util/time.cpp -L../external/lpj/lib -llpj -lbase -lsoil -ltree -lgrass -lnum -lclimate lpj.cpp -o lpj

# noah lsm
#gfortran -c --free-form --free-line-length-none ../external/noahlsm/module_model_constants.F
gfortran -g -c --free-form ~/桌面/simple_driver-v3.2/module_sf_noahlsm.F
gfortran -g -c --free-form ~/桌面/simple_driver-v3.2/module_Noahlsm_utility.F
gfortran -c --free-form --free-line-length-none ../external/noahlsm/kwm_date_utilities.F

g++ -I../include noahlsm.cpp ../src/model/noahlsm.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/model/model.cpp module_sf_noahlsm.o -o noahlsm/noahlsm -lgfortran
g++ -I../include noahlsm.cpp ../src/model/noahlsm.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/model/model.cpp /home/wlx/桌面/simple_driver-v3.2/*.o -o noahlsm/noahlsm -lgfortran

# colm with enkf
g++ -I./ -I../include colm_da.cpp ../src/model/commonlandmodel.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/util/landgrid.cpp ../src/model/model.cpp cldas.cpp -o colm/colmda -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/obs/radiancetransfermodel.cpp ../src/obs/qhmodel.cpp ../src/util/random.cpp ../src/da/enkf.cpp ../src/util/permittivity.cpp
# colm with enkf & openmp
g++ -I./ -I../include colm_da.cpp ../src/model/commonlandmodel.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/util/landgrid.cpp ../src/model/model.cpp cldas.cpp -o colm/colmda -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/obs/radiancetransfermodel.cpp ../src/obs/qhmodel.cpp ../src/util/random.cpp ../src/da/enkf.cpp ../src/util/permittivity.cpp -fopenmp
# colm with enkf & mpi
mpiCC -I./ -I../include colm_da_mpi.cpp ../src/model/commonlandmodel.cpp ../src/util/exception.cpp ../src/util/time.cpp ../src/util/landgrid.cpp ../src/model/model.cpp cldas_mpi.cpp -o colm/colmdampi -I../external/dSFMT  -DDSFMT_MEXP=19937 ../external/dSFMT/dSFMT.c ../src/obs/radiancetransfermodel.cpp ../src/obs/qhmodel.cpp ../src/util/random.cpp ../src/da/enkf.cpp ../src/util/permittivity.cpp

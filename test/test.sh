#time test
g++ -I../include timetest.cpp ../src/util/time.cpp -o timetest -lcpptest
g++ -I../include randomtest.cpp ../src/util/exception.cpp ../src/util/random.cpp ../external/dSFMT/dSFMT.c -lcpptest -I../external/dSFMT  -DDSFMT_MEXP=19937 -o randomtest
g++ -I../include permittivitytest.cpp ../src/util/exception.cpp ../src/util/permittivity.cpp -lcpptest -o permittivitytest
g++ -I../include gridtest.cpp ../src/util/exception.cpp ../src/util/landgrid.cpp -lcpptest -o gridtest

## aiem test
gfortran -c ../external/aiem/sigma.for
g++ -I../include aiemtest.cpp ../src/util/exception.cpp ../src/obs/radiancetransfermodel.cpp ../src/obs/aiem.cpp ../src/util/permittivity.cpp sigma.o -o aiemtest -lcpptest -lgfortran

## qhmodel test
g++ -I../include qhtest.cpp ../src/util/exception.cpp ../src/obs/radiancetransfermodel.cpp ../src/obs/qhmodel.cpp ../src/util/permittivity.cpp -o qhtest -lcpptest

## memls test
g++ -I../include memlstest.cpp ../src/util/exception.cpp ../src/obs/radiancetransfermodel.cpp ../src/obs/memls.cpp ../src/util/permittivity.cpp -o memlstest -lcpptest

## prosail test
rm *.o
gfortran -c ../external/prosail/*.f90
rm main_PROSAIL.o
g++ -lgfortran -I../include *.o prosailtest.cpp ../src/util/exception.cpp ../src/obs/radiancetransfermodel.cpp ../src/obs/prosail.cpp -o prosailtest -l cpptest

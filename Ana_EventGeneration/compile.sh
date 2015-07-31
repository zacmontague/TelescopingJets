
export PYTHIAPATH=$PWD/../../pythia8210

g++ MakeNTupleFromPythia.cc $PYTHIAPATH/lib/libpythia8.a -o MakeNTupleFromPythia -w -I$ROOTSYS/include -I$PYTHIAPATH/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath $PYTHIAPATH/lib -ldl -Wl,-rpath $ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs`

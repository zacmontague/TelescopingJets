export ROOTPATH=/Users/meehan/work/root
export PYTHIAPATH=/Users/meehan/work/TelescopingJets/pythia8210

g++ MakeNTupleFromPythia.cc $PYTHIAPATH/lib/libpythia8.a -o MakeNTupleFromPythia -w -I$ROOTPATH/include -I$PYTHIAPATH/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath $PYTHIAPATH/lib -ldl -Wl,-rpath $ROOTPATH/lib `$ROOTPATH/bin/root-config --glibs`

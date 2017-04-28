export PYTHIAPATH=/Users/meehan/work/TelescopingJets/TelescopingJets/pythia8223

echo
echo "Compiling with : "
echo "$ROOTSYS    : "${ROOTSYS}
echo "$PYTHIAPATH : "${PYTHIAPATH}
echo

g++ MakeNTupleFromPythia.cc $PYTHIAPATH/lib/libpythia8.a -o MakeNTupleFromPythia.exe -I$ROOTSYS/include  -I$PYTHIAPATH/include  -rpath $ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs` -std=c++11
export FJINSTALL=/Users/zacmon/Research/fastjet-3.3.0

echo
echo "Compiling with : "
echo "$ROOTSYS    : "${ROOTSYS}
echo "$FJINSTALL : "${FJINSTALL}
echo

g++ -o NTuplerTest.exe NTuplerTest.cc TelescopingJets.cc -I$ROOTSYS/include/root -I$FJINSTALL/include/fastjet  -rpath $ROOTSYS/lib `$ROOTSYS/bin/root-config --glibs` `$FJINSTALL/fastjet-config --cxxflags --libs --plugins` -lNsubjettiness -lEnergyCorrelator   -std=c++11

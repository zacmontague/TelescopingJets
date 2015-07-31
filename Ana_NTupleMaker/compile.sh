
export FJINSTALL=$PWD/../../fastjet-3.1.1

g++ -o NTupler NTupler.cc  TelescopingJets.cc -w -I$ROOTSYS/include -L$ROOTSYS/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad  -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,$ROOTSYS/lib -lm -ldl -pthread -m64 -I$ROOTSYS/include `$FJINSTALL/fastjet-config --cxxflags --libs --plugins` -lNsubjettiness
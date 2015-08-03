#compile
source compile.sh

#generate dijet
./MakeNTupleFromPythia 0 pythia_ttbar.root

#generate W'-->WZ
./MakeNTupleFromPythia 1 pythia_ttbar.root

#generate ttbar
./MakeNTupleFromPythia 2 pythia_ttbar.root
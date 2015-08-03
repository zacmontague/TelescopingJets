#compile
source compile.sh

#generate dijet
./MakeNTupleFromPythia 0 pythia_dijet.root

#generate W'-->WZ
./MakeNTupleFromPythia 1 pythia_wz.root

#generate ttbar
./MakeNTupleFromPythia 2 pythia_ttbar.root
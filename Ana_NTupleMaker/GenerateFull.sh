#compile
source compile.sh

#generate dijet
./NTupler 0 ../Ana_EventGeneration/pythia_dijet.root ntuple_dijet.root

#generate W'-->WZ
./NTupler 1 ../Ana_EventGeneration/pythia_wz.root ntuple_wz.root

#generate ttbar
./NTupler 2 ../Ana_EventGeneration/pythia_ttbar.root ntuple_ttbar.root
#Clean old root files
mkdir Data
rm Data/pythia_higgs_single.root
rm Data/ntuple_higgs_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 3 ../Data/pythia_higgs_single.root 1000
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 3 ../Data/pythia_higgs_single.root ../Data/ntuple_higgs_single.root
cd ..
#Clean old root files
mkdir Data
rm Data/pythia_ttbar_single.root
rm Data/ntuple_ttbar_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 1 ../Data/pythia_ttbar_single.root 1000
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 1 ../Data/pythia_ttbar_single.root ../Data/ntuple_ttbar_single.root
cd ..
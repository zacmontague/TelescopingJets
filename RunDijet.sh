#Clean old root files
mkdir Data
rm Data/pythia_dijet_single.root
rm Data/ntuple_dijet_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 0 ../Data/pythia_dijet_single.root 1000
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 0 ../Data/pythia_dijet_single.root ../Data/ntuple_dijet_single.root
cd ..
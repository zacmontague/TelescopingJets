#Clean old root files
mkdir Data
rm Data/pythia_wz_single.root
rm Data/ntuple_wz_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 2 ../Data/pythia_wz_single.root 1000
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 2 ../Data/pythia_wz_single.root ../Data/ntuple_wz_single.root
cd ..
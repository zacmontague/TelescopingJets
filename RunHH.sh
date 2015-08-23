#Clean old root files
mkdir Data
rm Data/pythia_hh_single.root
rm Data/ntuple_hh_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 4 ../Data/pythia_hh_single.root $1
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 4 ../Data/pythia_hh_single.root ../Data/ntuple_hh_single.root
cd ..
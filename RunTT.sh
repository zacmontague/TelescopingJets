#Clean old root files
mkdir Data
rm Data/pythia_tt_single.root
rm Data/ntuple_tt_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 3 ../Data/pythia_tt_single.root $1
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 3 ../Data/pythia_tt_single.root ../Data/ntuple_tt_single.root
cd ..
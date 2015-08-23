#Clean old root files
mkdir Data
rm Data/pythia_ww_single.root
rm Data/ntuple_ww_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 1 ../Data/pythia_ww_single.root $1
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 1 ../Data/pythia_ww_single.root ../Data/ntuple_ww_single.root
cd ..
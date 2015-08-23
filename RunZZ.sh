#Clean old root files
mkdir Data
rm Data/pythia_zz_single.root
rm Data/ntuple_zz_single.root

#Run generation
cd Ana_EventGeneration
./MakeNTupleFromPythia 2 ../Data/pythia_zz_single.root $1
cd ..

#Run ntupling
cd Ana_NTupleMaker
./NTupler 2 ../Data/pythia_zz_single.root ../Data/ntuple_zz_single.root
cd ..
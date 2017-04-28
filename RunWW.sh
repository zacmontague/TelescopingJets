#Clean old root files
mkdir Data
rm Data/pythia_ww.root
rm Data/ntuple_ww.root

nevents=1000

#Run generation
cd Ana_EventGeneration
source compile.sh
./MakeNTupleFromPythia.exe 1 ../Data/pythia_ww.root $nevents
cd ..

#Run ntupling
cd Ana_NTupleMaker
source compile.sh
./NTupler.exe 1 ../Data/pythia_ww.root ../Data/ntuple_ww.root
cd ..
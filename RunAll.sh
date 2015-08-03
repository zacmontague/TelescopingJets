#Clean old root files
mkdir Data
rm Data/*root
rm Ana_EventGeneration/*root
rm Ana_NTupleMaker/*root

#Run generation
cd Ana_EventGeneration
source GenerateFull.sh
cp *root ../Data/.
cd ..

#Run ntupling
cd Ana_NTupleMaker
source GenerateFull.sh
cp *root ../Data/.
cd ..

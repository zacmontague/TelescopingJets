echo "Running batch job with ..."
echo "Type:        "$1
echo "Name:        "$2
echo "NEvents:     "$3
echo "TransferDir: "$4

export TYPE=$1
export NAME=$2
export NEVENTS=$3
export TRANSFERDIR=$4

echo "Clean old root files"
mkdir Data
rm Data/pythia_$NAME.root
rm Data/ntuple_$NAME.root

echo "Run generation"
cd Ana_EventGeneration
./MakeNTupleFromPythia $TYPE ../Data/pythia_$NAME.root $NEVENTS
cd ..

echo "Run ntupling"
cd Ana_NTupleMaker
./NTupler $TYPE ../Data/pythia_$NAME.root ../Data/ntuple_$NAME.root
cd ..

echo "Transfer data back to where jobs were submitted"
echo "TransferDir: "$TRANSFERDIR
cp Data/pythia_$NAME.root $TRANSFERDIR/pythia_$NAME.root
cp Data/ntuple_$NAME.root $TRANSFERDIR/ntuple_$NAME.root
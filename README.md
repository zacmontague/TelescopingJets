# TelescopingJets
Analysis tools for the Telescoping Jets project

################################
# First time setup
################################
# After checking out the code which allows you to see this README, you must install three third party sources of code that are used by this outlined below with tips from our experience.  Note that these must be housed at the same level as the TelescopingJets code that you checked out.  This is particularly true in the case of event generation with Pythia8 

# Pythia8 : http://home.thep.lu.se/~torbjorn/Pythia.html
1) No tips here.  But if you are not familiar with Pythia then it is worthwhile working through the worksheet outlined here http://home.thep.lu.se/~torbjorn/pythia8/worksheet8183.pdf since the EventGeneration code is built starting from this.

# fastjet : http://fastjet.fr/
1) Check out the latest stable version - I have it working with fastjet-3.1.1 and the paths are configured to look for this, so if you don't check out this version, you will need to modify some paths in the compile.sh scripts

# fastjet contrib : http://fastjet.hepforge.org/contrib/
1) When installing, it is necessary to point the code to the fastjet install directory as shown in the install directions.  This is done via the [--fastjet-config=FILE] argument.  For me, on Mac OSX, when I installed fastjet, this meant pointing it to the fastjet-config at ../fastjet-3.1.1/fastjet-config.  If you do this correctly, then you should see the /contrib directory housed in /usr/local/include/fastjet/.


################################
# Running EventGeneration
################################
Where : Ana_EventGeneration
How : 
$ source compile.sh
$ ./MakeNTupleFromPythia 1 pythia_wz.root

This will generate 500 events and save to the pythia_wz.root output file a flat ntuple that contains truth vectors of the relevant particles (W,Z) for truth tagging later along with vector<double> branches for all of the truth particles with final state particles.

################################
# Running NTupling
################################
Where : Ana_NTupleMaker
How : 
$ source compile.sh
$ ./NTupler 1 pythia\_wz.root NTuple\_wz.root

This will process the 500 events generated with the MakeNTupleFromPythia code and do the following:
1) run jet finding on the truth constituents 
2) determine truth flavor of jets via dR matching using the truth 4-vectors from the (W,Z)
3) calculate basic TJet volatility for these jets

The output will be saved in the file NTuple_wz.root

################################
# Running MiniNTupleAnalysis
################################
Where : Ana_MiniNTupleAnalysis
How : 
$ python MakeROCS.py

As it currently stands, this requires that you obtain an additional directory from outside of GitHub that contains the root files that will be read in.  This directory is housed on DropBox at:

https://www.dropbox.com/sh/hfu09io74gswafp/AACf_8KucJmvTsDB-1iQ_JeDa?dl=0

and you should copy this entire "Data" directory to the same level as the repository you just checked out.  



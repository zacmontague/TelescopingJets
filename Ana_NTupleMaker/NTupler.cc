/*----------------------------------------------------------------------

TITLE : NTupler.cc

DESCRIPTION : Takes input final state particle ntuple from Ana_EventGeneration
and outputs a flat ntuple that contains jet properties created via fastjet
in this code to be used for post-analysis using Ana_MiniNTupleAnalysis. (NOTE!!)
All of the TelescopingJets code from the fastjet contrib should be contained within this code

COMPILE :
$ source compile.sh

RUN :
$ ./NTupler <type> <input> <output>

type   : 0 = dijet , 1 = W'->WZ , 2 = ttbar
input  : Input file from Ana_MiniNTupleAnalysis
output : Anything you want - but being logical

//----------------------------------------------------------------------*/


#include "NTupler.h"

int main(int argc, char* argv[]){


  //exit if you dont pass a run card
  if(argc<4){
    cout<<"You need to specify more arguments"<<endl;
    return 1;
  }

  //process
  int ProcessType = atoi(argv[1]);

  //inputfile
  string InputFile = argv[2];

  //outputfile
  string OutputFile = argv[3];

  //debug flag
  bool debug=false;

  //dR truth matching
  double dR_match = 1.0;

  //////////////////////////////////////////////
  //INPUT
  //////////////////////////////////////////////
  //get input file and tree
  filein = new TFile( InputFile.c_str() );
  treein = (TTree*)filein->Get( "tree" );
  treein->Print();

  //set up branch linking to addresses
  treein->SetBranchAddress("fspart_id", &fspart_id);
  treein->SetBranchAddress("fspart_pt", &fspart_pt);
  treein->SetBranchAddress("fspart_eta",&fspart_eta);
  treein->SetBranchAddress("fspart_phi",&fspart_phi);
  treein->SetBranchAddress("fspart_m",  &fspart_m);

  treein->SetBranchAddress("truth_q1_pt",  &truth_q1_pt);
  treein->SetBranchAddress("truth_q1_eta", &truth_q1_eta);
  treein->SetBranchAddress("truth_q1_phi", &truth_q1_phi);
  treein->SetBranchAddress("truth_q1_m",   &truth_q1_m);

  treein->SetBranchAddress("truth_q2_pt",  &truth_q2_pt);
  treein->SetBranchAddress("truth_q2_eta", &truth_q2_eta);
  treein->SetBranchAddress("truth_q2_phi", &truth_q2_phi);
  treein->SetBranchAddress("truth_q2_m",   &truth_q2_m);

  treein->SetBranchAddress("truth_t1_pt",  &truth_t1_pt);
  treein->SetBranchAddress("truth_t1_eta", &truth_t1_eta);
  treein->SetBranchAddress("truth_t1_phi", &truth_t1_phi);
  treein->SetBranchAddress("truth_t1_m",   &truth_t1_m);

  treein->SetBranchAddress("truth_t2_pt",  &truth_t2_pt);
  treein->SetBranchAddress("truth_t2_eta", &truth_t2_eta);
  treein->SetBranchAddress("truth_t2_phi", &truth_t2_phi);
  treein->SetBranchAddress("truth_t2_m",   &truth_t2_m);

  treein->SetBranchAddress("truth_W_pt",  &truth_W_pt);
  treein->SetBranchAddress("truth_W_eta", &truth_W_eta);
  treein->SetBranchAddress("truth_W_phi", &truth_W_phi);
  treein->SetBranchAddress("truth_W_m",   &truth_W_m);

  treein->SetBranchAddress("truth_Z_pt",  &truth_Z_pt);
  treein->SetBranchAddress("truth_Z_eta", &truth_Z_eta);
  treein->SetBranchAddress("truth_Z_phi", &truth_Z_phi);
  treein->SetBranchAddress("truth_Z_m",   &truth_Z_m);


  //////////////////////////////////////////////
  //OUTPUT
  //////////////////////////////////////////////
  fileout = new TFile( OutputFile.c_str() ,"RECREATE");

  treeout = new TTree("JetTree","JetTree");

  treeout->Branch("NumberOfVertices",    &NumberOfVertices);

  treeout->Branch("TruthRaw_flavor",        &TruthRaw_flavor);
  treeout->Branch("TruthRaw_pt",            &TruthRaw_pt);
  treeout->Branch("TruthRaw_eta",           &TruthRaw_eta);
  treeout->Branch("TruthRaw_phi",           &TruthRaw_phi);
  treeout->Branch("TruthRaw_m",             &TruthRaw_m);
  treeout->Branch("TruthRaw_Tau1",          &TruthRaw_Tau1);
  treeout->Branch("TruthRaw_Tau2",          &TruthRaw_Tau2);
  treeout->Branch("TruthRaw_Tau21",         &TruthRaw_Tau21);
  treeout->Branch("TruthRaw_C2",            &TruthRaw_C2);
  treeout->Branch("TruthRaw_D2",            &TruthRaw_D2);
  treeout->Branch("TruthRaw_C3",            &TruthRaw_C3);
  treeout->Branch("TruthRaw_TJet_m1",       &TruthRaw_TJet_m1);
  treeout->Branch("TruthRaw_TJet_m2",       &TruthRaw_TJet_m2);
  treeout->Branch("TruthRaw_TJet_Tau1",     &TruthRaw_TJet_Tau1);
  treeout->Branch("TruthRaw_TJet_Tau2",     &TruthRaw_TJet_Tau2);
  treeout->Branch("TruthRaw_TJet_Tau21",    &TruthRaw_TJet_Tau21);
  treeout->Branch("TruthRaw_TJet_C2",       &TruthRaw_TJet_C2);
  treeout->Branch("TruthRaw_TJet_D2",       &TruthRaw_TJet_D2);
  treeout->Branch("TruthRaw_TJet_C3",       &TruthRaw_TJet_C3);

  treeout->Branch("TruthPileup_flavor",        &TruthPileup_flavor);
  treeout->Branch("TruthPileup_pt",            &TruthPileup_pt);
  treeout->Branch("TruthPileup_eta",           &TruthPileup_eta);
  treeout->Branch("TruthPileup_phi",           &TruthPileup_phi);
  treeout->Branch("TruthPileup_m",             &TruthPileup_m);
  treeout->Branch("TruthPileup_Tau1",          &TruthPileup_Tau1);
  treeout->Branch("TruthPileup_Tau2",          &TruthPileup_Tau2);
  treeout->Branch("TruthPileup_Tau21",         &TruthPileup_Tau21);
  treeout->Branch("TruthPileup_C2",            &TruthPileup_C2);
  treeout->Branch("TruthPileup_D2",            &TruthPileup_D2);
  treeout->Branch("TruthPileup_C3",            &TruthPileup_C3);
  treeout->Branch("TruthPileup_TJet_m1",       &TruthPileup_TJet_m1);
  treeout->Branch("TruthPileup_TJet_m2",       &TruthPileup_TJet_m2);
  treeout->Branch("TruthPileup_TJet_Tau1",     &TruthPileup_TJet_Tau1);
  treeout->Branch("TruthPileup_TJet_Tau2",     &TruthPileup_TJet_Tau2);
  treeout->Branch("TruthPileup_TJet_Tau21",    &TruthPileup_TJet_Tau21);
  treeout->Branch("TruthPileup_TJet_C2",       &TruthPileup_TJet_C2);
  treeout->Branch("TruthPileup_TJet_D2",       &TruthPileup_TJet_D2);
  treeout->Branch("TruthPileup_TJet_C3",       &TruthPileup_TJet_C3);

  treeout->Branch("RecoRaw_flavor",        &RecoRaw_flavor);
  treeout->Branch("RecoRaw_pt",            &RecoRaw_pt);
  treeout->Branch("RecoRaw_eta",           &RecoRaw_eta);
  treeout->Branch("RecoRaw_phi",           &RecoRaw_phi);
  treeout->Branch("RecoRaw_m",             &RecoRaw_m);
  treeout->Branch("RecoRaw_Tau1",          &RecoRaw_Tau1);
  treeout->Branch("RecoRaw_Tau2",          &RecoRaw_Tau2);
  treeout->Branch("RecoRaw_Tau21",         &RecoRaw_Tau21);
  treeout->Branch("RecoRaw_C2",            &RecoRaw_C2);
  treeout->Branch("RecoRaw_D2",            &RecoRaw_D2);
  treeout->Branch("RecoRaw_C3",            &RecoRaw_C3);
  treeout->Branch("RecoRaw_TJet_m1",       &RecoRaw_TJet_m1);
  treeout->Branch("RecoRaw_TJet_m2",       &RecoRaw_TJet_m2);
  treeout->Branch("RecoRaw_TJet_Tau1",     &RecoRaw_TJet_Tau1);
  treeout->Branch("RecoRaw_TJet_Tau2",     &RecoRaw_TJet_Tau2);
  treeout->Branch("RecoRaw_TJet_Tau21",    &RecoRaw_TJet_Tau21);
  treeout->Branch("RecoRaw_TJet_C2",       &RecoRaw_TJet_C2);
  treeout->Branch("RecoRaw_TJet_D2",       &RecoRaw_TJet_D2);
  treeout->Branch("RecoRaw_TJet_C3",       &RecoRaw_TJet_C3);

  treeout->Branch("RecoPileup_flavor",        &RecoPileup_flavor);
  treeout->Branch("RecoPileup_pt",            &RecoPileup_pt);
  treeout->Branch("RecoPileup_eta",           &RecoPileup_eta);
  treeout->Branch("RecoPileup_phi",           &RecoPileup_phi);
  treeout->Branch("RecoPileup_m",             &RecoPileup_m);
  treeout->Branch("RecoPileup_Tau1",          &RecoPileup_Tau1);
  treeout->Branch("RecoPileup_Tau2",          &RecoPileup_Tau2);
  treeout->Branch("RecoPileup_Tau21",         &RecoPileup_Tau21);
  treeout->Branch("RecoPileup_C2",            &RecoPileup_C2);
  treeout->Branch("RecoPileup_D2",            &RecoPileup_D2);
  treeout->Branch("RecoPileup_C3",            &RecoPileup_C3);
  treeout->Branch("RecoPileup_TJet_m1",       &RecoPileup_TJet_m1);
  treeout->Branch("RecoPileup_TJet_m2",       &RecoPileup_TJet_m2);
  treeout->Branch("RecoPileup_TJet_Tau1",     &RecoPileup_TJet_Tau1);
  treeout->Branch("RecoPileup_TJet_Tau2",     &RecoPileup_TJet_Tau2);
  treeout->Branch("RecoPileup_TJet_Tau21",    &RecoPileup_TJet_Tau21);
  treeout->Branch("RecoPileup_TJet_C2",       &RecoPileup_TJet_C2);
  treeout->Branch("RecoPileup_TJet_D2",       &RecoPileup_TJet_D2);
  treeout->Branch("RecoPileup_TJet_C3",       &RecoPileup_TJet_C3);

  treeout->Branch("TruthRawTrim_flavor",        &TruthRawTrim_flavor);
  treeout->Branch("TruthRawTrim_pt",            &TruthRawTrim_pt);
  treeout->Branch("TruthRawTrim_eta",           &TruthRawTrim_eta);
  treeout->Branch("TruthRawTrim_phi",           &TruthRawTrim_phi);
  treeout->Branch("TruthRawTrim_m",             &TruthRawTrim_m);
  treeout->Branch("TruthRawTrim_Tau1",          &TruthRawTrim_Tau1);
  treeout->Branch("TruthRawTrim_Tau2",          &TruthRawTrim_Tau2);
  treeout->Branch("TruthRawTrim_Tau21",         &TruthRawTrim_Tau21);
  treeout->Branch("TruthRawTrim_C2",            &TruthRawTrim_C2);
  treeout->Branch("TruthRawTrim_D2",            &TruthRawTrim_D2);
  treeout->Branch("TruthRawTrim_C3",            &TruthRawTrim_C3);
  treeout->Branch("TruthRawTrim_TJet_m1",       &TruthRawTrim_TJet_m1);
  treeout->Branch("TruthRawTrim_TJet_m2",       &TruthRawTrim_TJet_m2);
  treeout->Branch("TruthRawTrim_TJet_Tau1",     &TruthRawTrim_TJet_Tau1);
  treeout->Branch("TruthRawTrim_TJet_Tau2",     &TruthRawTrim_TJet_Tau2);
  treeout->Branch("TruthRawTrim_TJet_Tau21",    &TruthRawTrim_TJet_Tau21);
  treeout->Branch("TruthRawTrim_TJet_C2",       &TruthRawTrim_TJet_C2);
  treeout->Branch("TruthRawTrim_TJet_D2",       &TruthRawTrim_TJet_D2);
  treeout->Branch("TruthRawTrim_TJet_C3",       &TruthRawTrim_TJet_C3);

  treeout->Branch("TruthPileupTrim_flavor",        &TruthPileupTrim_flavor);
  treeout->Branch("TruthPileupTrim_pt",            &TruthPileupTrim_pt);
  treeout->Branch("TruthPileupTrim_eta",           &TruthPileupTrim_eta);
  treeout->Branch("TruthPileupTrim_phi",           &TruthPileupTrim_phi);
  treeout->Branch("TruthPileupTrim_m",             &TruthPileupTrim_m);
  treeout->Branch("TruthPileupTrim_Tau1",          &TruthPileupTrim_Tau1);
  treeout->Branch("TruthPileupTrim_Tau2",          &TruthPileupTrim_Tau2);
  treeout->Branch("TruthPileupTrim_Tau21",         &TruthPileupTrim_Tau21);
  treeout->Branch("TruthPileupTrim_C2",            &TruthPileupTrim_C2);
  treeout->Branch("TruthPileupTrim_D2",            &TruthPileupTrim_D2);
  treeout->Branch("TruthPileupTrim_C3",            &TruthPileupTrim_C3);
  treeout->Branch("TruthPileupTrim_TJet_m1",       &TruthPileupTrim_TJet_m1);
  treeout->Branch("TruthPileupTrim_TJet_m2",       &TruthPileupTrim_TJet_m2);
  treeout->Branch("TruthPileupTrim_TJet_Tau1",     &TruthPileupTrim_TJet_Tau1);
  treeout->Branch("TruthPileupTrim_TJet_Tau2",     &TruthPileupTrim_TJet_Tau2);
  treeout->Branch("TruthPileupTrim_TJet_Tau21",    &TruthPileupTrim_TJet_Tau21);
  treeout->Branch("TruthPileupTrim_TJet_C2",       &TruthPileupTrim_TJet_C2);
  treeout->Branch("TruthPileupTrim_TJet_D2",       &TruthPileupTrim_TJet_D2);
  treeout->Branch("TruthPileupTrim_TJet_C3",       &TruthPileupTrim_TJet_C3);

  treeout->Branch("RecoRawTrim_flavor",        &RecoRawTrim_flavor);
  treeout->Branch("RecoRawTrim_pt",            &RecoRawTrim_pt);
  treeout->Branch("RecoRawTrim_eta",           &RecoRawTrim_eta);
  treeout->Branch("RecoRawTrim_phi",           &RecoRawTrim_phi);
  treeout->Branch("RecoRawTrim_m",             &RecoRawTrim_m);
  treeout->Branch("RecoRawTrim_Tau1",          &RecoRawTrim_Tau1);
  treeout->Branch("RecoRawTrim_Tau2",          &RecoRawTrim_Tau2);
  treeout->Branch("RecoRawTrim_Tau21",         &RecoRawTrim_Tau21);
  treeout->Branch("RecoRawTrim_C2",            &RecoRawTrim_C2);
  treeout->Branch("RecoRawTrim_D2",            &RecoRawTrim_D2);
  treeout->Branch("RecoRawTrim_C3",            &RecoRawTrim_C3);
  treeout->Branch("RecoRawTrim_TJet_m1",       &RecoRawTrim_TJet_m1);
  treeout->Branch("RecoRawTrim_TJet_m2",       &RecoRawTrim_TJet_m2);
  treeout->Branch("RecoRawTrim_TJet_Tau1",     &RecoRawTrim_TJet_Tau1);
  treeout->Branch("RecoRawTrim_TJet_Tau2",     &RecoRawTrim_TJet_Tau2);
  treeout->Branch("RecoRawTrim_TJet_Tau21",    &RecoRawTrim_TJet_Tau21);
  treeout->Branch("RecoRawTrim_TJet_C2",       &RecoRawTrim_TJet_C2);
  treeout->Branch("RecoRawTrim_TJet_D2",       &RecoRawTrim_TJet_D2);
  treeout->Branch("RecoRawTrim_TJet_C3",       &RecoRawTrim_TJet_C3);

  treeout->Branch("RecoPileupTrim_flavor",        &RecoPileupTrim_flavor);
  treeout->Branch("RecoPileupTrim_pt",            &RecoPileupTrim_pt);
  treeout->Branch("RecoPileupTrim_eta",           &RecoPileupTrim_eta);
  treeout->Branch("RecoPileupTrim_phi",           &RecoPileupTrim_phi);
  treeout->Branch("RecoPileupTrim_m",             &RecoPileupTrim_m);
  treeout->Branch("RecoPileupTrim_Tau1",          &RecoPileupTrim_Tau1);
  treeout->Branch("RecoPileupTrim_Tau2",          &RecoPileupTrim_Tau2);
  treeout->Branch("RecoPileupTrim_Tau21",         &RecoPileupTrim_Tau21);
  treeout->Branch("RecoPileupTrim_C2",            &RecoPileupTrim_C2);
  treeout->Branch("RecoPileupTrim_D2",            &RecoPileupTrim_D2);
  treeout->Branch("RecoPileupTrim_C3",            &RecoPileupTrim_C3);
  treeout->Branch("RecoPileupTrim_TJet_m1",       &RecoPileupTrim_TJet_m1);
  treeout->Branch("RecoPileupTrim_TJet_m2",       &RecoPileupTrim_TJet_m2);
  treeout->Branch("RecoPileupTrim_TJet_Tau1",     &RecoPileupTrim_TJet_Tau1);
  treeout->Branch("RecoPileupTrim_TJet_Tau2",     &RecoPileupTrim_TJet_Tau2);
  treeout->Branch("RecoPileupTrim_TJet_Tau21",    &RecoPileupTrim_TJet_Tau21);
  treeout->Branch("RecoPileupTrim_TJet_C2",       &RecoPileupTrim_TJet_C2);
  treeout->Branch("RecoPileupTrim_TJet_D2",       &RecoPileupTrim_TJet_D2);
  treeout->Branch("RecoPileupTrim_TJet_C3",       &RecoPileupTrim_TJet_C3);

  ////////////////////////////////
  //random number generator for pileup
  ////////////////////////////////
  TRandom3 *rand_pileup = new TRandom3();


  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  cout<<"Number of events: "<<nEvents<<endl;

  for (Long64_t jentry=0; jentry<nEvents; jentry++) {

    if(jentry%500==0)
      cout<<"Event: "<<jentry<<endl;

    //Get next event from input ntuple
    filein->cd();
    treein->GetEntry(jentry);

    /////////////////////////////
    //Reset branches for next event
    /////////////////////////////
    ResetBranches();

    ///////////////////////////////////////////////////
    //read in all final state particles for jet building from pythia input
    ///////////////////////////////////////////////////
    vector<PseudoJet> input_particles;
    input_particles.clear();

    int n_fspart = fspart_id->size();
    for(int i_fspart=0; i_fspart<n_fspart; i_fspart++){

      if(debug){
        cout<<fspart_id->at(i_fspart)<<"  "
            <<fspart_pt->at(i_fspart)<<"  "
            <<fspart_eta->at(i_fspart)<<"  "
            <<fspart_phi->at(i_fspart)<<"  "
            <<fspart_m->at(i_fspart)<<"  "<<endl;
      }

      TLorentzVector temp_p4;
      temp_p4.SetPtEtaPhiM(fspart_pt ->at(i_fspart),
                           fspart_eta->at(i_fspart),
                           fspart_phi->at(i_fspart),
                           fspart_m  ->at(i_fspart));


      input_particles.push_back(PseudoJet(temp_p4.Px(),temp_p4.Py(),temp_p4.Pz(),temp_p4.E()));
    }

    //////////////////////////////////////////////
    //make new input particles collection with pileup
    //////////////////////////////////////////////
    vector<PseudoJet> input_particles_Pileup;
    input_particles_Pileup.clear();
    for(int ipart=0; ipart<n_fspart; ipart++){
      input_particles_Pileup.push_back(input_particles.at(ipart));
    }

    int n_pileup_vertices      = (int)rand_pileup->Poisson(10);
    int n_particles_per_vertex = 5;
    int n_pileup_particles = n_pileup_vertices*n_particles_per_vertex;

    NumberOfVertices = n_pileup_vertices;

    if(debug) cout<<"Pileup: "<<NumberOfVertices<<"  "<<n_particles_per_vertex<<"  "<<n_pileup_particles<<endl;

    for(int ipart=0; ipart<n_pileup_particles; ipart++){

      double m  = 0.0;
      double px = rand_pileup->Gaus(0,5.0);
      double py = rand_pileup->Gaus(0,5.0);
      double pz = rand_pileup->Gaus(0,5.0);
      double E  = pow( m*m + px*px + py*py + pz*pz , 0.5);

      if(debug) cout<<"Pileup: "<<ipart<<"  "<<px<<"  "<<py<<"  "<<pz<<"  "<<E<<endl;

      input_particles_Pileup.push_back(PseudoJet(px,py,pz,E));

    }


    //////////////////////////////////////////////
    //make pseudocalorimeter cells
    //////////////////////////////////////////////
    vector<PseudoJet> calo_cells        = ToyCalorimeter(input_particles);
    vector<PseudoJet> calo_cells_Pileup = ToyCalorimeter(input_particles_Pileup);


    //////////////////////////////////////////////
    // get the resulting jets ordered in pt
    //////////////////////////////////////////////
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);

    fastjet::ClusterSequence clust_seq_TruthRaw(input_particles, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_TruthPileup(input_particles_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthPileup = sorted_by_pt(clust_seq_TruthPileup.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoRaw(calo_cells, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoRaw = sorted_by_pt(clust_seq_RecoRaw.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoPileup(calo_cells_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoPileup = sorted_by_pt(clust_seq_RecoPileup.inclusive_jets(5.0));



    if(debug){
      // label the columns
      cout<<"jet#  pt  eta  phi  mass"<<endl;
      cout<<"Inclusive"<<endl;
      // print out the details for each jet
      for (unsigned int i = 0; i < inclusive_jets_TruthRaw.size(); i++) {
        cout<<i<<"  "<<inclusive_jets_TruthRaw[i].pt()
               <<"  "<<inclusive_jets_TruthRaw[i].eta()
               <<"  "<<inclusive_jets_TruthRaw[i].phi()
               <<"  "<<inclusive_jets_TruthRaw[i].m()<<endl;
      }
    }

    //////////////////////////////////////////////
    //Setup tools for substructure calculation
    //////////////////////////////////////////////

    //Telescoping jets
    fastjet::contrib::KT_Axes axes_def;
    std::vector<double> r_values;
    for(int i=0; i < 20; i++){
      r_values.push_back( 0.1+i*(0.6-0.1)/(20-1) );
    }
    TelescopingJets T_Mass(axes_def,r_values);

    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);

    //Energy correlation functions
    fastjet::contrib::EnergyCorrelatorC2 ecfC2(1.);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1.);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1.);

    // Filtering with a pt cut as for trimming (arXiv:0912.1342)
    Transformer *trimmer = new Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.05) );
    const Transformer &f = *trimmer;

    /////////////////////////////////////////////
    //Get truth objects for truth matching
    /////////////////////////////////////////////
    TLorentzVector truth_q1;
    truth_q1.SetPtEtaPhiM(truth_q1_pt,truth_q1_eta,truth_q1_phi,truth_q1_m);
    TLorentzVector truth_q2;
    truth_q2.SetPtEtaPhiM(truth_q2_pt,truth_q2_eta,truth_q2_phi,truth_q2_m);
    TLorentzVector truth_W;
    truth_W.SetPtEtaPhiM(truth_W_pt,truth_W_eta,truth_W_phi,truth_W_m);
    TLorentzVector truth_Z;
    truth_Z.SetPtEtaPhiM(truth_Z_pt,truth_Z_eta,truth_Z_phi,truth_Z_m);
    TLorentzVector truth_t1;
    truth_t1.SetPtEtaPhiM(truth_t1_pt,truth_t1_eta,truth_t1_phi,truth_t1_m);
    TLorentzVector truth_t2;
    truth_t2.SetPtEtaPhiM(truth_t2_pt,truth_t2_eta,truth_t2_phi,truth_t2_m);




    /////////////////////////////
    //TruthRaw
    /////////////////////////////
    if(debug) cout<<"TruthRaw jet"<<endl;
    for(int ijet=0; ijet<inclusive_jets_TruthRaw.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets_TruthRaw.at(ijet).pt(),
                           inclusive_jets_TruthRaw.at(ijet).eta(),
                           inclusive_jets_TruthRaw.at(ijet).phi(),
                           inclusive_jets_TruthRaw.at(ijet).m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(inclusive_jets_TruthRaw[ijet]);
      tempJet_Tau2           = nsub2(inclusive_jets_TruthRaw[ijet]);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(inclusive_jets_TruthRaw[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_TruthRaw[ijet]);
      tempJet_C3             = ecfC3(inclusive_jets_TruthRaw[ijet]);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_TruthRaw[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_TruthRaw[ijet]);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, inclusive_jets_TruthRaw[ijet], 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, inclusive_jets_TruthRaw[ijet], 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, inclusive_jets_TruthRaw[ijet], 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(inclusive_jets_TruthRaw[ijet], 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(inclusive_jets_TruthRaw[ijet], 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(inclusive_jets_TruthRaw[ijet], 0.1, 2.);

      TruthRaw_flavor    .push_back(tempJet_flavor);
      TruthRaw_pt        .push_back(tempJet_pt);
      TruthRaw_eta       .push_back(tempJet_eta);
      TruthRaw_phi       .push_back(tempJet_phi);
      TruthRaw_m         .push_back(tempJet_m);
      TruthRaw_Tau1      .push_back(tempJet_Tau1);
      TruthRaw_Tau2      .push_back(tempJet_Tau2);
      TruthRaw_Tau21     .push_back(tempJet_Tau21);
      TruthRaw_C2        .push_back(tempJet_C2);
      TruthRaw_D2        .push_back(tempJet_D2);
      TruthRaw_C3        .push_back(tempJet_C3);
      TruthRaw_TJet_m1   .push_back(tempJet_TJet_m1);
      TruthRaw_TJet_m2   .push_back(tempJet_TJet_m2);
      TruthRaw_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      TruthRaw_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      TruthRaw_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      TruthRaw_TJet_C2   .push_back(tempJet_TJet_C2);
      TruthRaw_TJet_D2   .push_back(tempJet_TJet_D2);
      TruthRaw_TJet_C3   .push_back(tempJet_TJet_C3);
    }


    /////////////////////////////
    //TruthPileup
    /////////////////////////////
    if(debug) cout<<"TruthPileup jet"<<endl;
    for(int ijet=0; ijet<inclusive_jets_TruthPileup.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets_TruthPileup.at(ijet).pt(),
                           inclusive_jets_TruthPileup.at(ijet).eta(),
                           inclusive_jets_TruthPileup.at(ijet).phi(),
                           inclusive_jets_TruthPileup.at(ijet).m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(inclusive_jets_TruthPileup[ijet]);
      tempJet_Tau2           = nsub2(inclusive_jets_TruthPileup[ijet]);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(inclusive_jets_TruthPileup[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_TruthPileup[ijet]);
      tempJet_C3             = ecfC3(inclusive_jets_TruthPileup[ijet]);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_TruthPileup[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_TruthPileup[ijet]);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, inclusive_jets_TruthPileup[ijet], 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, inclusive_jets_TruthPileup[ijet], 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, inclusive_jets_TruthPileup[ijet], 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(inclusive_jets_TruthPileup[ijet], 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(inclusive_jets_TruthPileup[ijet], 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(inclusive_jets_TruthPileup[ijet], 0.1, 2.);

      TruthPileup_flavor    .push_back(tempJet_flavor);
      TruthPileup_pt        .push_back(tempJet_pt);
      TruthPileup_eta       .push_back(tempJet_eta);
      TruthPileup_phi       .push_back(tempJet_phi);
      TruthPileup_m         .push_back(tempJet_m);
      TruthPileup_Tau1      .push_back(tempJet_Tau1);
      TruthPileup_Tau2      .push_back(tempJet_Tau2);
      TruthPileup_Tau21     .push_back(tempJet_Tau21);
      TruthPileup_C2        .push_back(tempJet_C2);
      TruthPileup_D2        .push_back(tempJet_D2);
      TruthPileup_C3        .push_back(tempJet_C3);
      TruthPileup_TJet_m1   .push_back(tempJet_TJet_m1);
      TruthPileup_TJet_m2   .push_back(tempJet_TJet_m2);
      TruthPileup_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      TruthPileup_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      TruthPileup_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      TruthPileup_TJet_C2   .push_back(tempJet_TJet_C2);
      TruthPileup_TJet_D2   .push_back(tempJet_TJet_D2);
      TruthPileup_TJet_C3   .push_back(tempJet_TJet_C3);

    }


    /////////////////////////////
    //RecoRaw
    /////////////////////////////
    if(debug) cout<<"RecoRaw jet"<<endl;
    for(int ijet=0; ijet<inclusive_jets_RecoRaw.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets_RecoRaw.at(ijet).pt(),
                           inclusive_jets_RecoRaw.at(ijet).eta(),
                           inclusive_jets_RecoRaw.at(ijet).phi(),
                           inclusive_jets_RecoRaw.at(ijet).m());

      /////////////////////////////////
      //Getting Reco label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(inclusive_jets_RecoRaw[ijet]);
      tempJet_Tau2           = nsub2(inclusive_jets_RecoRaw[ijet]);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(inclusive_jets_RecoRaw[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_RecoRaw[ijet]);
      tempJet_C3             = ecfC3(inclusive_jets_RecoRaw[ijet]);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_RecoRaw[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_RecoRaw[ijet]);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, inclusive_jets_RecoRaw[ijet], 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, inclusive_jets_RecoRaw[ijet], 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, inclusive_jets_RecoRaw[ijet], 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(inclusive_jets_RecoRaw[ijet], 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(inclusive_jets_RecoRaw[ijet], 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(inclusive_jets_RecoRaw[ijet], 0.1, 2.);

      RecoRaw_flavor    .push_back(tempJet_flavor);
      RecoRaw_pt        .push_back(tempJet_pt);
      RecoRaw_eta       .push_back(tempJet_eta);
      RecoRaw_phi       .push_back(tempJet_phi);
      RecoRaw_m         .push_back(tempJet_m);
      RecoRaw_Tau1      .push_back(tempJet_Tau1);
      RecoRaw_Tau2      .push_back(tempJet_Tau2);
      RecoRaw_Tau21     .push_back(tempJet_Tau21);
      RecoRaw_C2        .push_back(tempJet_C2);
      RecoRaw_D2        .push_back(tempJet_D2);
      RecoRaw_C3        .push_back(tempJet_C3);
      RecoRaw_TJet_m1   .push_back(tempJet_TJet_m1);
      RecoRaw_TJet_m2   .push_back(tempJet_TJet_m2);
      RecoRaw_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      RecoRaw_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      RecoRaw_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      RecoRaw_TJet_C2   .push_back(tempJet_TJet_C2);
      RecoRaw_TJet_D2   .push_back(tempJet_TJet_D2);
      RecoRaw_TJet_C3   .push_back(tempJet_TJet_C3);

    }


    /////////////////////////////
    //RecoPileup
    /////////////////////////////
    if(debug) cout<<"RecoPileup jet"<<endl;
    for(int ijet=0; ijet<inclusive_jets_RecoPileup.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets_RecoPileup.at(ijet).pt(),
                           inclusive_jets_RecoPileup.at(ijet).eta(),
                           inclusive_jets_RecoPileup.at(ijet).phi(),
                           inclusive_jets_RecoPileup.at(ijet).m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(inclusive_jets_RecoPileup[ijet]);
      tempJet_Tau2           = nsub2(inclusive_jets_RecoPileup[ijet]);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(inclusive_jets_RecoPileup[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_RecoPileup[ijet]);
      tempJet_C3             = ecfC3(inclusive_jets_RecoPileup[ijet]);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_RecoPileup[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_RecoPileup[ijet]);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, inclusive_jets_RecoPileup[ijet], 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, inclusive_jets_RecoPileup[ijet], 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, inclusive_jets_RecoPileup[ijet], 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(inclusive_jets_RecoPileup[ijet], 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(inclusive_jets_RecoPileup[ijet], 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(inclusive_jets_RecoPileup[ijet], 0.1, 2.);

      RecoPileup_flavor    .push_back(tempJet_flavor);
      RecoPileup_pt        .push_back(tempJet_pt);
      RecoPileup_eta       .push_back(tempJet_eta);
      RecoPileup_phi       .push_back(tempJet_phi);
      RecoPileup_m         .push_back(tempJet_m);
      RecoPileup_Tau1      .push_back(tempJet_Tau1);
      RecoPileup_Tau2      .push_back(tempJet_Tau2);
      RecoPileup_Tau21     .push_back(tempJet_Tau21);
      RecoPileup_C2        .push_back(tempJet_C2);
      RecoPileup_D2        .push_back(tempJet_D2);
      RecoPileup_C3        .push_back(tempJet_C3);
      RecoPileup_TJet_m1   .push_back(tempJet_TJet_m1);
      RecoPileup_TJet_m2   .push_back(tempJet_TJet_m2);
      RecoPileup_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      RecoPileup_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      RecoPileup_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      RecoPileup_TJet_C2   .push_back(tempJet_TJet_C2);
      RecoPileup_TJet_D2   .push_back(tempJet_TJet_D2);
      RecoPileup_TJet_C3   .push_back(tempJet_TJet_C3);

    }



    /////////////////////////////
    //TruthRawTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_TruthRaw.size(); ijet++) {
      PseudoJet groomed_jet = f(inclusive_jets_TruthRaw[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;


      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(groomed_jet);
      tempJet_Tau2           = nsub2(groomed_jet);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);
      tempJet_C3             = ecfC3(groomed_jet);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, groomed_jet, 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, groomed_jet, 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, groomed_jet, 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(groomed_jet, 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(groomed_jet, 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(groomed_jet, 0.1, 2.);

      TruthRawTrim_flavor    .push_back(tempJet_flavor);
      TruthRawTrim_pt        .push_back(tempJet_pt);
      TruthRawTrim_eta       .push_back(tempJet_eta);
      TruthRawTrim_phi       .push_back(tempJet_phi);
      TruthRawTrim_m         .push_back(tempJet_m);
      TruthRawTrim_Tau1      .push_back(tempJet_Tau1);
      TruthRawTrim_Tau2      .push_back(tempJet_Tau2);
      TruthRawTrim_Tau21     .push_back(tempJet_Tau21);
      TruthRawTrim_C2        .push_back(tempJet_C2);
      TruthRawTrim_D2        .push_back(tempJet_D2);
      TruthRawTrim_C3        .push_back(tempJet_C3);
      TruthRawTrim_TJet_m1   .push_back(tempJet_TJet_m1);
      TruthRawTrim_TJet_m2   .push_back(tempJet_TJet_m2);
      TruthRawTrim_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      TruthRawTrim_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      TruthRawTrim_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      TruthRawTrim_TJet_C2   .push_back(tempJet_TJet_C2);
      TruthRawTrim_TJet_D2   .push_back(tempJet_TJet_D2);
      TruthRawTrim_TJet_C3   .push_back(tempJet_TJet_C3);

    }


    /////////////////////////////
    //TruthPileupTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_TruthPileup.size(); ijet++) {
      PseudoJet groomed_jet = f(inclusive_jets_TruthPileup[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;


      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(groomed_jet);
      tempJet_Tau2           = nsub2(groomed_jet);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);
      tempJet_C3             = ecfC3(groomed_jet);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, groomed_jet, 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, groomed_jet, 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, groomed_jet, 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(groomed_jet, 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(groomed_jet, 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(groomed_jet, 0.1, 2.);

      TruthPileupTrim_flavor    .push_back(tempJet_flavor);
      TruthPileupTrim_pt        .push_back(tempJet_pt);
      TruthPileupTrim_eta       .push_back(tempJet_eta);
      TruthPileupTrim_phi       .push_back(tempJet_phi);
      TruthPileupTrim_m         .push_back(tempJet_m);
      TruthPileupTrim_Tau1      .push_back(tempJet_Tau1);
      TruthPileupTrim_Tau2      .push_back(tempJet_Tau2);
      TruthPileupTrim_Tau21     .push_back(tempJet_Tau21);
      TruthPileupTrim_C2        .push_back(tempJet_C2);
      TruthPileupTrim_D2        .push_back(tempJet_D2);
      TruthPileupTrim_C3        .push_back(tempJet_C3);
      TruthPileupTrim_TJet_m1   .push_back(tempJet_TJet_m1);
      TruthPileupTrim_TJet_m2   .push_back(tempJet_TJet_m2);
      TruthPileupTrim_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      TruthPileupTrim_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      TruthPileupTrim_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      TruthPileupTrim_TJet_C2   .push_back(tempJet_TJet_C2);
      TruthPileupTrim_TJet_D2   .push_back(tempJet_TJet_D2);
      TruthPileupTrim_TJet_C3   .push_back(tempJet_TJet_C3);

    }

    /////////////////////////////
    //RecoRawTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_RecoRaw.size(); ijet++) {
      PseudoJet groomed_jet = f(inclusive_jets_RecoRaw[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;


      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(groomed_jet);
      tempJet_Tau2           = nsub2(groomed_jet);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);
      tempJet_C3             = ecfC3(groomed_jet);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, groomed_jet, 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, groomed_jet, 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, groomed_jet, 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(groomed_jet, 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(groomed_jet, 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(groomed_jet, 0.1, 2.);

      RecoRawTrim_flavor    .push_back(tempJet_flavor);
      RecoRawTrim_pt        .push_back(tempJet_pt);
      RecoRawTrim_eta       .push_back(tempJet_eta);
      RecoRawTrim_phi       .push_back(tempJet_phi);
      RecoRawTrim_m         .push_back(tempJet_m);
      RecoRawTrim_Tau1      .push_back(tempJet_Tau1);
      RecoRawTrim_Tau2      .push_back(tempJet_Tau2);
      RecoRawTrim_Tau21     .push_back(tempJet_Tau21);
      RecoRawTrim_C2        .push_back(tempJet_C2);
      RecoRawTrim_D2        .push_back(tempJet_D2);
      RecoRawTrim_C3        .push_back(tempJet_C3);
      RecoRawTrim_TJet_m1   .push_back(tempJet_TJet_m1);
      RecoRawTrim_TJet_m2   .push_back(tempJet_TJet_m2);
      RecoRawTrim_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      RecoRawTrim_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      RecoRawTrim_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      RecoRawTrim_TJet_C2   .push_back(tempJet_TJet_C2);
      RecoRawTrim_TJet_D2   .push_back(tempJet_TJet_D2);
      RecoRawTrim_TJet_C3   .push_back(tempJet_TJet_C3);

    }

    /////////////////////////////
    //RecoPileupTrim
    /////////////////////////////
    for (unsigned int ijet = 0; ijet < inclusive_jets_RecoPileup.size(); ijet++) {
      PseudoJet groomed_jet = f(inclusive_jets_RecoPileup[ijet]);

      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(groomed_jet.pt(),
                           groomed_jet.eta(),
                           groomed_jet.phi(),
                           groomed_jet.m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      if(debug){
        cout<<"DeltaR: "<<endl
            <<"q1:  "<<jettemp.DeltaR(truth_q1)<<endl
            <<"q2:  "<<jettemp.DeltaR(truth_q2)<<endl
            <<"W:   "<<jettemp.DeltaR(truth_W)<<endl
            <<"Z:   "<<jettemp.DeltaR(truth_Z)<<endl
            <<"t1:  "<<jettemp.DeltaR(truth_t1)<<endl
            <<"t2:  "<<jettemp.DeltaR(truth_t2)<<endl;
      }
      int jetflavor = -1;
      if(jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){ jetflavor = 0; }
      else if(jettemp.DeltaR(truth_W)<dR_match){  jetflavor = 1; }
      else if(jettemp.DeltaR(truth_Z)<dR_match){  jetflavor = 2; }
      else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){ jetflavor = 3; }
      else{ jetflavor = -1; }

      if(debug) cout<<"FillingJet: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;


      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau1           = nsub1(groomed_jet);
      tempJet_Tau2           = nsub2(groomed_jet);
      if(tempJet_Tau1!=0)
        tempJet_Tau21        = tempJet_Tau2/tempJet_Tau1;
      tempJet_C2             = ecfC2(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);
      tempJet_C3             = ecfC3(groomed_jet);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);
      tempJet_TJet_Tau1      = T_Nsubjettiness(1, groomed_jet, 1., 2.);
      tempJet_TJet_Tau2      = T_Nsubjettiness(2, groomed_jet, 1., 2.);
      tempJet_TJet_Tau21     = T_NsubjettinessRatio(2, 1, groomed_jet, 1., 2.);
      tempJet_TJet_C2        = T_EnergyCorrelator_C2(groomed_jet, 0.1, 2.);
      tempJet_TJet_D2        = T_EnergyCorrelator_D2(groomed_jet, 0.1, 2.);
      tempJet_TJet_C3        = T_EnergyCorrelator_C3(groomed_jet, 0.1, 2.);

      RecoPileupTrim_flavor    .push_back(tempJet_flavor);
      RecoPileupTrim_pt        .push_back(tempJet_pt);
      RecoPileupTrim_eta       .push_back(tempJet_eta);
      RecoPileupTrim_phi       .push_back(tempJet_phi);
      RecoPileupTrim_m         .push_back(tempJet_m);
      RecoPileupTrim_Tau1      .push_back(tempJet_Tau1);
      RecoPileupTrim_Tau2      .push_back(tempJet_Tau2);
      RecoPileupTrim_Tau21     .push_back(tempJet_Tau21);
      RecoPileupTrim_C2        .push_back(tempJet_C2);
      RecoPileupTrim_D2        .push_back(tempJet_D2);
      RecoPileupTrim_C3        .push_back(tempJet_C3);
      RecoPileupTrim_TJet_m1   .push_back(tempJet_TJet_m1);
      RecoPileupTrim_TJet_m2   .push_back(tempJet_TJet_m2);
      RecoPileupTrim_TJet_Tau1 .push_back(tempJet_TJet_Tau1);
      RecoPileupTrim_TJet_Tau2 .push_back(tempJet_TJet_Tau2);
      RecoPileupTrim_TJet_Tau21.push_back(tempJet_TJet_Tau21);
      RecoPileupTrim_TJet_C2   .push_back(tempJet_TJet_C2);
      RecoPileupTrim_TJet_D2   .push_back(tempJet_TJet_D2);
      RecoPileupTrim_TJet_C3   .push_back(tempJet_TJet_C3);

    }









    //////////////////////////////////////
    //Fill event into tree
    //////////////////////////////////////
    if(debug) cout<<"Filling Tree"<<endl;
    treeout->Fill();
  }


  /////////////////////////////////
  //Write the output TTree to the OutputFile
  /////////////////////////////////
  fileout->cd();
  treeout->Write();
  fileout->Close();

  return 0;

}



///=========================================
/// Reset Branches
///=========================================
void ResetBranches(){

  NumberOfVertices = 0;

  TruthRaw_flavor.clear();
  TruthRaw_pt.clear();
  TruthRaw_eta.clear();
  TruthRaw_phi.clear();
  TruthRaw_m.clear();
  TruthRaw_Tau1.clear();
  TruthRaw_Tau2.clear();
  TruthRaw_Tau21.clear();
  TruthRaw_C2.clear();
  TruthRaw_D2.clear();
  TruthRaw_C3.clear();
  TruthRaw_TJet_m1.clear();
  TruthRaw_TJet_m2.clear();
  TruthRaw_TJet_Tau1.clear();
  TruthRaw_TJet_Tau2.clear();
  TruthRaw_TJet_Tau21.clear();
  TruthRaw_TJet_C2.clear();
  TruthRaw_TJet_D2.clear();
  TruthRaw_TJet_C3.clear();

  TruthPileup_flavor.clear();
  TruthPileup_pt.clear();
  TruthPileup_eta.clear();
  TruthPileup_phi.clear();
  TruthPileup_m.clear();
  TruthPileup_Tau1.clear();
  TruthPileup_Tau2.clear();
  TruthPileup_Tau21.clear();
  TruthPileup_C2.clear();
  TruthPileup_D2.clear();
  TruthPileup_C3.clear();
  TruthPileup_TJet_m1.clear();
  TruthPileup_TJet_m2.clear();
  TruthPileup_TJet_Tau1.clear();
  TruthPileup_TJet_Tau2.clear();
  TruthPileup_TJet_Tau21.clear();
  TruthPileup_TJet_C2.clear();
  TruthPileup_TJet_D2.clear();
  TruthPileup_TJet_C3.clear();

  RecoRaw_flavor.clear();
  RecoRaw_pt.clear();
  RecoRaw_eta.clear();
  RecoRaw_phi.clear();
  RecoRaw_m.clear();
  RecoRaw_Tau1.clear();
  RecoRaw_Tau2.clear();
  RecoRaw_Tau21.clear();
  RecoRaw_C2.clear();
  RecoRaw_D2.clear();
  RecoRaw_C3.clear();
  RecoRaw_TJet_m1.clear();
  RecoRaw_TJet_m2.clear();
  RecoRaw_TJet_Tau1.clear();
  RecoRaw_TJet_Tau2.clear();
  RecoRaw_TJet_Tau21.clear();
  RecoRaw_TJet_C2.clear();
  RecoRaw_TJet_D2.clear();
  RecoRaw_TJet_C3.clear();

  RecoPileup_flavor.clear();
  RecoPileup_pt.clear();
  RecoPileup_eta.clear();
  RecoPileup_phi.clear();
  RecoPileup_m.clear();
  RecoPileup_Tau1.clear();
  RecoPileup_Tau2.clear();
  RecoPileup_Tau21.clear();
  RecoPileup_C2.clear();
  RecoPileup_D2.clear();
  RecoPileup_C3.clear();
  RecoPileup_TJet_m1.clear();
  RecoPileup_TJet_m2.clear();
  RecoPileup_TJet_Tau1.clear();
  RecoPileup_TJet_Tau2.clear();
  RecoPileup_TJet_Tau21.clear();
  RecoPileup_TJet_C2.clear();
  RecoPileup_TJet_D2.clear();
  RecoPileup_TJet_C3.clear();



  TruthRawTrim_flavor.clear();
  TruthRawTrim_pt.clear();
  TruthRawTrim_eta.clear();
  TruthRawTrim_phi.clear();
  TruthRawTrim_m.clear();
  TruthRawTrim_Tau1.clear();
  TruthRawTrim_Tau2.clear();
  TruthRawTrim_Tau21.clear();
  TruthRawTrim_C2.clear();
  TruthRawTrim_D2.clear();
  TruthRawTrim_C3.clear();
  TruthRawTrim_TJet_m1.clear();
  TruthRawTrim_TJet_m2.clear();
  TruthRawTrim_TJet_Tau1.clear();
  TruthRawTrim_TJet_Tau2.clear();
  TruthRawTrim_TJet_Tau21.clear();
  TruthRawTrim_TJet_C2.clear();
  TruthRawTrim_TJet_D2.clear();
  TruthRawTrim_TJet_C3.clear();

  TruthPileupTrim_flavor.clear();
  TruthPileupTrim_pt.clear();
  TruthPileupTrim_eta.clear();
  TruthPileupTrim_phi.clear();
  TruthPileupTrim_m.clear();
  TruthPileupTrim_Tau1.clear();
  TruthPileupTrim_Tau2.clear();
  TruthPileupTrim_Tau21.clear();
  TruthPileupTrim_C2.clear();
  TruthPileupTrim_D2.clear();
  TruthPileupTrim_C3.clear();
  TruthPileupTrim_TJet_m1.clear();
  TruthPileupTrim_TJet_m2.clear();
  TruthPileupTrim_TJet_Tau1.clear();
  TruthPileupTrim_TJet_Tau2.clear();
  TruthPileupTrim_TJet_Tau21.clear();
  TruthPileupTrim_TJet_C2.clear();
  TruthPileupTrim_TJet_D2.clear();
  TruthPileupTrim_TJet_C3.clear();

  RecoRawTrim_flavor.clear();
  RecoRawTrim_pt.clear();
  RecoRawTrim_eta.clear();
  RecoRawTrim_phi.clear();
  RecoRawTrim_m.clear();
  RecoRawTrim_Tau1.clear();
  RecoRawTrim_Tau2.clear();
  RecoRawTrim_Tau21.clear();
  RecoRawTrim_C2.clear();
  RecoRawTrim_D2.clear();
  RecoRawTrim_C3.clear();
  RecoRawTrim_TJet_m1.clear();
  RecoRawTrim_TJet_m2.clear();
  RecoRawTrim_TJet_Tau1.clear();
  RecoRawTrim_TJet_Tau2.clear();
  RecoRawTrim_TJet_Tau21.clear();
  RecoRawTrim_TJet_C2.clear();
  RecoRawTrim_TJet_D2.clear();
  RecoRawTrim_TJet_C3.clear();

  RecoPileupTrim_flavor.clear();
  RecoPileupTrim_pt.clear();
  RecoPileupTrim_eta.clear();
  RecoPileupTrim_phi.clear();
  RecoPileupTrim_m.clear();
  RecoPileupTrim_Tau1.clear();
  RecoPileupTrim_Tau2.clear();
  RecoPileupTrim_Tau21.clear();
  RecoPileupTrim_C2.clear();
  RecoPileupTrim_D2.clear();
  RecoPileupTrim_C3.clear();
  RecoPileupTrim_TJet_m1.clear();
  RecoPileupTrim_TJet_m2.clear();
  RecoPileupTrim_TJet_Tau1.clear();
  RecoPileupTrim_TJet_Tau2.clear();
  RecoPileupTrim_TJet_Tau21.clear();
  RecoPileupTrim_TJet_C2.clear();
  RecoPileupTrim_TJet_D2.clear();
  RecoPileupTrim_TJet_C3.clear();

}


///=========================================
/// Calorimeter Simulation
///=========================================
vector<PseudoJet> ToyCalorimeter(vector<PseudoJet> truth_particles) {
  const double pi = 3.14159265359;
  const double etaLim = 5.0;
  const int nEta = 100;
  const int nPhi = 63;
  double dEta = 2*etaLim/nEta;
  double dPhi = 2*pi/nPhi;

  double tower[nEta][nPhi];
  for (int i = 0; i < nEta; i++)  for (int j = 0; j < nPhi; j++)  tower[i][j] = -0.001;

  vector<fastjet::PseudoJet> cell_particles;
  for (int p=0; p < (int)truth_particles.size(); p++) {
    fastjet::PseudoJet part = truth_particles.at(p);

    int etaCell = int((part.eta()+etaLim)/dEta);
    int phiCell = int(part.phi()/dPhi);
    if (etaCell >= 0 && etaCell < nEta && phiCell >=0 && phiCell < nPhi){
      tower[etaCell][phiCell] += part.e();
    }
  }

  for (int i = 0; i < nEta; i++)  for (int j = 0; j < nPhi; j++) {
    if (tower[i][j] > 0) {
      double etaLocal = -etaLim + dEta*(i+0.5);
      double phiLocal = dPhi*(j+0.5);
      double thetaLocal = 2*atan(exp(-etaLocal));
      cell_particles.push_back(fastjet::PseudoJet(sin(thetaLocal)*cos(phiLocal),sin(thetaLocal)*sin(phiLocal),cos(thetaLocal),1)*tower[i][j]);
    }
  }
  return cell_particles;
}


///=========================================
/// Telescoping N-subjettiness
///=========================================
double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max) {
  vector<double> taus; taus.clear();
  for (int i = 0; i < 20; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(20-1);
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    taus.push_back(nsub(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(taus);
}

double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max) {
  vector<double> taus; taus.clear();
  for (int i = 0; i < 20; i++) {

    double beta = beta_min + i*(beta_max - beta_min)/(20-1);

    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);

    fastjet::contrib::Nsubjettiness nsub_num(N_num, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub_den(N_den, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);

    double num=nsub_num(input);
    double den=nsub_den(input);

    if(den!=0)
      taus.push_back(num/den);
    else
      taus.push_back(-1.0);

  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(taus);
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < 20; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(20-1);
    fastjet::contrib::EnergyCorrelatorC2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
}

double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < 20; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(20-1);
    fastjet::contrib::EnergyCorrelatorD2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
}

double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < 20; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(20-1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  return getVolatility(ecfs);
}
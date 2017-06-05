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

type   : 0 = dijet , 1 = G*->W+W- , 2 = ttbar
input  : Input file from Ana_MiniNTupleAnalysis
output : Anything you want - but being logical

//----------------------------------------------------------------------*/


#include "NTupler.h"

int main(int argc, char* argv[]){


  //exit if you dont pass a run card
  if(argc<4){
    cout<<"You need to specify more arguments"<<endl;
    cout<<"Arg1 = process type"<<endl;
    cout<<"Arg2 = input file path and name"<<endl;
    cout<<"Arg3 = output file path and name"<<endl;
    cout<<"Arg4 = debug flag (optional)"<<endl;
    return 1;
  }

  //process
  string ProcessType = argv[1];

  //inputfile
  string InputFile = argv[2];

  //outputfile
  string OutputFile = argv[3];

  //debug flag
  bool debug=false;
  if(argc>=5){
    string argdebug = argv[4];
    if(argdebug=="debug")
      debug=true;
  }

  //print out the input arguments
  cout<<"InputArguments:  ProcessType="<<ProcessType<<"   InputFile="<<InputFile<<"  OutputFile="<<OutputFile<<"  Debug="<<debug<<endl;

  //dR truth matching
  dR_match = 1.0;

  //////////////////////////////////////////////
  //INPUT
  //////////////////////////////////////////////
  //get input file and tree
  filein = new TFile( InputFile.c_str() );
  treein = (TTree*)filein->Get( "tree" );
  if(debug) treein->Print();

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

  treein->SetBranchAddress("truth_W1_pt",  &truth_W1_pt);
  treein->SetBranchAddress("truth_W1_eta", &truth_W1_eta);
  treein->SetBranchAddress("truth_W1_phi", &truth_W1_phi);
  treein->SetBranchAddress("truth_W1_m",   &truth_W1_m);

  treein->SetBranchAddress("truth_W2_pt",  &truth_W2_pt);
  treein->SetBranchAddress("truth_W2_eta", &truth_W2_eta);
  treein->SetBranchAddress("truth_W2_phi", &truth_W2_phi);
  treein->SetBranchAddress("truth_W2_m",   &truth_W2_m);

  treein->SetBranchAddress("truth_H_pt",  &truth_H_pt);
  treein->SetBranchAddress("truth_H_eta", &truth_H_eta);
  treein->SetBranchAddress("truth_H_phi", &truth_H_phi);
  treein->SetBranchAddress("truth_H_m",   &truth_H_m);

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
  treeout->Branch("TruthRaw_Tau21",         &TruthRaw_Tau21);
  treeout->Branch("TruthRaw_Tau32",         &TruthRaw_Tau32);
  treeout->Branch("TruthRaw_D2",            &TruthRaw_D2);
  treeout->Branch("TruthRaw_T2jet_angle",   &TruthRaw_T2jet_angle);
  treeout->Branch("TruthRaw_T2jet",         &TruthRaw_T2jet);
  treeout->Branch("TruthRaw_T3jet_angle",   &TruthRaw_T3jet_angle);
  treeout->Branch("TruthRaw_T3jet",         &TruthRaw_T3jet);
  treeout->Branch("TruthRaw_Tpruning",      &TruthRaw_Tpruning);
  treeout->Branch("TruthRaw_Ttrimming",     &TruthRaw_Ttrimming);
  treeout->Branch("TruthRaw_Taktreclustering",	&TruthRaw_Taktreclustering);
  treeout->Branch("TruthRaw_Tktreclustering",	&TruthRaw_Tktreclustering);
  treeout->Branch("TruthRaw_TJet_m1",       &TruthRaw_TJet_m1);
  treeout->Branch("TruthRaw_TJet_m2",       &TruthRaw_TJet_m2);

  treeout->Branch("TruthRawTrim_flavor",        &TruthRawTrim_flavor);
  treeout->Branch("TruthRawTrim_pt",            &TruthRawTrim_pt);
  treeout->Branch("TruthRawTrim_eta",           &TruthRawTrim_eta);
  treeout->Branch("TruthRawTrim_phi",           &TruthRawTrim_phi);
  treeout->Branch("TruthRawTrim_m",             &TruthRawTrim_m);
  treeout->Branch("TruthRawTrim_Tau21",         &TruthRawTrim_Tau21);
  treeout->Branch("TruthRawTrim_Tau32",         &TruthRawTrim_Tau32);
  treeout->Branch("TruthRawTrim_D2",            &TruthRawTrim_D2);
  treeout->Branch("TruthRawTrim_T2jet_angle",   &TruthRawTrim_T2jet_angle);
  treeout->Branch("TruthRawTrim_T2jet",         &TruthRawTrim_T2jet);
  treeout->Branch("TruthRawTrim_T3jet_angle",   &TruthRawTrim_T3jet_angle);
  treeout->Branch("TruthRawTrim_T3jet",         &TruthRawTrim_T3jet);
  treeout->Branch("TruthRawTrim_Tpruning",      &TruthRawTrim_Tpruning);
  treeout->Branch("TruthRawTrim_Ttrimming",     &TruthRawTrim_Ttrimming);
  treeout->Branch("TruthRawTrim_Taktreclustering",	&TruthRawTrim_Taktreclustering);
  treeout->Branch("TruthRawTrim_Tktreclustering",	&TruthRawTrim_Tktreclustering);
  treeout->Branch("TruthRawTrim_TJet_m1",       &TruthRawTrim_TJet_m1);
  treeout->Branch("TruthRawTrim_TJet_m2",       &TruthRawTrim_TJet_m2);

  treeout->Branch("RecoRaw_flavor",        &RecoRaw_flavor);
  treeout->Branch("RecoRaw_pt",            &RecoRaw_pt);
  treeout->Branch("RecoRaw_eta",           &RecoRaw_eta);
  treeout->Branch("RecoRaw_phi",           &RecoRaw_phi);
  treeout->Branch("RecoRaw_m",             &RecoRaw_m);
  treeout->Branch("RecoRaw_Tau21",         &RecoRaw_Tau21);
  treeout->Branch("RecoRaw_Tau32",         &RecoRaw_Tau32);
  treeout->Branch("RecoRaw_D2",            &RecoRaw_D2);
  treeout->Branch("RecoRaw_T2jet_angle",   &RecoRaw_T2jet_angle);
  treeout->Branch("RecoRaw_T2jet",         &RecoRaw_T2jet);
  treeout->Branch("RecoRaw_T3jet_angle",   &RecoRaw_T3jet_angle);
  treeout->Branch("RecoRaw_T3jet",         &RecoRaw_T3jet);
  treeout->Branch("RecoRaw_Tpruning",      &RecoRaw_Tpruning);
  treeout->Branch("RecoRaw_Ttrimming",     &RecoRaw_Ttrimming);
  treeout->Branch("RecoRaw_Taktreclustering",	&RecoRaw_Taktreclustering);
  treeout->Branch("RecoRaw_Tktreclustering",	&RecoRaw_Tktreclustering);
  treeout->Branch("RecoRaw_TJet_m1",       &RecoRaw_TJet_m1);
  treeout->Branch("RecoRaw_TJet_m2",       &RecoRaw_TJet_m2);

  treeout->Branch("RecoRawTrim_flavor",        &RecoRawTrim_flavor);
  treeout->Branch("RecoRawTrim_pt",            &RecoRawTrim_pt);
  treeout->Branch("RecoRawTrim_eta",           &RecoRawTrim_eta);
  treeout->Branch("RecoRawTrim_phi",           &RecoRawTrim_phi);
  treeout->Branch("RecoRawTrim_m",             &RecoRawTrim_m);
  treeout->Branch("RecoRawTrim_Tau21",         &RecoRawTrim_Tau21);
  treeout->Branch("RecoRawTrim_Tau32",         &RecoRawTrim_Tau32);
  treeout->Branch("RecoRawTrim_D2",            &RecoRawTrim_D2);
  treeout->Branch("RecoRawTrim_T2jet_angle",   &RecoRawTrim_T2jet_angle);
  treeout->Branch("RecoRawTrim_T2jet",         &RecoRawTrim_T2jet);
  treeout->Branch("RecoRawTrim_T3jet_angle",   &RecoRawTrim_T3jet_angle);
  treeout->Branch("RecoRawTrim_T3jet",         &RecoRawTrim_T3jet);
  treeout->Branch("RecoRawTrim_Tpruning",      &RecoRawTrim_Tpruning);
  treeout->Branch("RecoRawTrim_Ttrimming",     &RecoRawTrim_Ttrimming);
  treeout->Branch("RecoRawTrim_Taktreclustering",	&RecoRawTrim_Taktreclustering);
  treeout->Branch("RecoRawTrim_Tktreclustering",	&RecoRawTrim_Tktreclustering);
  treeout->Branch("RecoRawTrim_TJet_m1",       &RecoRawTrim_TJet_m1);
  treeout->Branch("RecoRawTrim_TJet_m2",       &RecoRawTrim_TJet_m2);

  //////////////////////////////////////////////
  //random number generator for pileup
  //////////////////////////////////////////////
  TRandom3 *rand_pileup = new TRandom3();


  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  cout<<"Number of events: "<<nEvents<<endl;

  for (Long64_t jentry=0; jentry<nEvents; jentry++) {

    if(jentry%10==0)
      cout<<"NTupler: ProcessType="<<ProcessType<<"  entry="<<jentry<<endl;

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

/*
    //////////////////////////////////////////////
    //make new input particles collection with pileup
    //////////////////////////////////////////////
    //this will be using min bias events from simulations

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
*/

    //////////////////////////////////////////////
    //make pseudocalorimeter cells
    //////////////////////////////////////////////
    vector<PseudoJet> calo_cells        = ToyCalorimeter(input_particles);
    //vector<PseudoJet> calo_cells_Pileup = ToyCalorimeter(input_particles_Pileup);


    //////////////////////////////////////////////
    // get the resulting jets ordered in pt
    //////////////////////////////////////////////
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.0);

    fastjet::ClusterSequence clust_seq_TruthRaw(input_particles, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthRaw = sorted_by_pt(clust_seq_TruthRaw.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoRaw(calo_cells, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoRaw = sorted_by_pt(clust_seq_RecoRaw.inclusive_jets(5.0));
/*
    fastjet::ClusterSequence clust_seq_TruthPileup(input_particles_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_TruthPileup = sorted_by_pt(clust_seq_TruthPileup.inclusive_jets(5.0));

    fastjet::ClusterSequence clust_seq_RecoPileup(calo_cells_Pileup, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_RecoPileup = sorted_by_pt(clust_seq_RecoPileup.inclusive_jets(5.0));
*/


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

    //Telescoping jets (this looks like the Telescoping reclustering)
    fastjet::contrib::KT_Axes axes_def;
    std::vector<double> r_values;
    int N_r = 20;
    double r_min = 0.1;
    double r_max = 0.6;
    for(int i=0; i < N_r; i++){
      r_values.push_back( r_min+i*(r_max-r_min)/(N_r-1) );
    }
    TelescopingJets T_Mass(axes_def,r_values);

    //Energy correlation functions
    fastjet::contrib::EnergyCorrelatorC2 ecfC2(1.);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1.);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC3(2, 1.);

    // Filtering with a pt cut as for trimming (arXiv:0912.1342)
    double Rfilt0 = 0.3;
    double fcut0 = 0.05;
    Transformer *trimmer = new Filter(JetDefinition(kt_algorithm, Rfilt0), SelectorPtFractionMin(fcut0) );
    const Transformer &f = *trimmer;

    /////////////////////////////////////////////
    //Get truth objects for truth matching
    /////////////////////////////////////////////
    truth_q1.SetPtEtaPhiM(truth_q1_pt,truth_q1_eta,truth_q1_phi,truth_q1_m);
    truth_q2.SetPtEtaPhiM(truth_q2_pt,truth_q2_eta,truth_q2_phi,truth_q2_m);
    truth_t1.SetPtEtaPhiM(truth_t1_pt,truth_t1_eta,truth_t1_phi,truth_t1_m);
    truth_t2.SetPtEtaPhiM(truth_t2_pt,truth_t2_eta,truth_t2_phi,truth_t2_m);
    truth_W1.SetPtEtaPhiM(truth_W1_pt,truth_W1_eta,truth_W1_phi,truth_W1_m);
    truth_W2.SetPtEtaPhiM(truth_W2_pt,truth_W2_eta,truth_W2_phi,truth_W2_m);
    truth_H.SetPtEtaPhiM(truth_H_pt,truth_H_eta,truth_H_phi,truth_H_m);

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
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Raw   : flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      if(jetflavor==-1)
        continue;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(inclusive_jets_TruthRaw[ijet]);
      tempJet_Tau32          = GetTau32(inclusive_jets_TruthRaw[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_TruthRaw[ijet]);

      TSub  T2SubOutput     = T_2Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutput.min_angle;
      tempJet_T2jet          = T2SubOutput.volatility;

      TSub  T3SubOutput     = T_3Subjet(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutput.min_angle;
      tempJet_T3jet          = T3SubOutput.volatility;

      tempJet_Tpruning       = T_Pruning (inclusive_jets_TruthRaw[ijet], 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(inclusive_jets_TruthRaw[ijet], 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(inclusive_jets_TruthRaw[ijet], 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_TruthRaw[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_TruthRaw[ijet]);

      TruthRaw_flavor     .push_back(tempJet_flavor);
      TruthRaw_pt         .push_back(tempJet_pt);
      TruthRaw_eta        .push_back(tempJet_eta);
      TruthRaw_phi        .push_back(tempJet_phi);
      TruthRaw_m          .push_back(tempJet_m);
      TruthRaw_Tau21      .push_back(tempJet_Tau21);
      TruthRaw_Tau32      .push_back(tempJet_Tau32);
      TruthRaw_D2         .push_back(tempJet_D2);
      TruthRaw_T2jet_angle.push_back(tempJet_T2jet_angle);
      TruthRaw_T2jet      .push_back(tempJet_T2jet);
      TruthRaw_T3jet_angle.push_back(tempJet_T3jet_angle);
      TruthRaw_T3jet      .push_back(tempJet_T3jet);
      TruthRaw_Tpruning   .push_back(tempJet_Tpruning);
      TruthRaw_Ttrimming  .push_back(tempJet_Ttrimming);
      TruthRaw_Taktreclustering.push_back(tempJet_Taktreclustering);
      TruthRaw_Tktreclustering.push_back(tempJet_Tktreclustering);
      TruthRaw_TJet_m1    .push_back(tempJet_TJet_m1);
      TruthRaw_TJet_m2    .push_back(tempJet_TJet_m2);
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
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Trimmed: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      if(jetflavor==-1)
        continue;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(groomed_jet);
      tempJet_Tau32          = GetTau32(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);

      TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutputTrim.min_angle;
      tempJet_T2jet          = T2SubOutputTrim.volatility;

      TSub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutputTrim.min_angle;
      tempJet_T3jet          = T3SubOutputTrim.volatility;

      tempJet_Tpruning       = T_Pruning (groomed_jet, 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(groomed_jet, 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);

      TruthRawTrim_flavor     .push_back(tempJet_flavor);
      TruthRawTrim_pt         .push_back(tempJet_pt);
      TruthRawTrim_eta        .push_back(tempJet_eta);
      TruthRawTrim_phi        .push_back(tempJet_phi);
      TruthRawTrim_m          .push_back(tempJet_m);
      TruthRawTrim_Tau21      .push_back(tempJet_Tau21);
      TruthRawTrim_Tau32      .push_back(tempJet_Tau32);
      TruthRawTrim_D2         .push_back(tempJet_D2);
      TruthRawTrim_T2jet_angle.push_back(tempJet_T2jet_angle);
      TruthRawTrim_T2jet      .push_back(tempJet_T2jet);
      TruthRawTrim_T3jet_angle.push_back(tempJet_T3jet_angle);
      TruthRawTrim_T3jet      .push_back(tempJet_T3jet);
      TruthRawTrim_Tpruning   .push_back(tempJet_Tpruning);
      TruthRawTrim_Ttrimming  .push_back(tempJet_Ttrimming);
      TruthRawTrim_Taktreclustering .push_back(tempJet_Taktreclustering);
      TruthRawTrim_Tktreclustering .push_back(tempJet_Tktreclustering);
      TruthRawTrim_TJet_m1    .push_back(tempJet_TJet_m1);
      TruthRawTrim_TJet_m2    .push_back(tempJet_TJet_m2);
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
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Raw   : flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      if(jetflavor==-1)
        continue;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(inclusive_jets_RecoRaw[ijet]);
      tempJet_Tau32          = GetTau32(inclusive_jets_RecoRaw[ijet]);
      tempJet_D2             = ecfD2(inclusive_jets_RecoRaw[ijet]);

      TSub  T2SubOutput     = T_2Subjet(inclusive_jets_RecoRaw[ijet], 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutput.min_angle;
      tempJet_T2jet          = T2SubOutput.volatility;

      TSub  T3SubOutput     = T_3Subjet(inclusive_jets_RecoRaw[ijet], 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutput.min_angle;
      tempJet_T3jet          = T3SubOutput.volatility;

      tempJet_Tpruning       = T_Pruning (inclusive_jets_RecoRaw[ijet], 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(inclusive_jets_RecoRaw[ijet], 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(inclusive_jets_RecoRaw[ijet], 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(inclusive_jets_RecoRaw[ijet], 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,inclusive_jets_RecoRaw[ijet]);
      tempJet_TJet_m2        = T_Mass(2,inclusive_jets_RecoRaw[ijet]);

      RecoRaw_flavor     .push_back(tempJet_flavor);
      RecoRaw_pt         .push_back(tempJet_pt);
      RecoRaw_eta        .push_back(tempJet_eta);
      RecoRaw_phi        .push_back(tempJet_phi);
      RecoRaw_m          .push_back(tempJet_m);
      RecoRaw_Tau21      .push_back(tempJet_Tau21);
      RecoRaw_Tau32      .push_back(tempJet_Tau32);
      RecoRaw_D2         .push_back(tempJet_D2);
      RecoRaw_T2jet_angle.push_back(tempJet_T2jet_angle);
      RecoRaw_T2jet      .push_back(tempJet_T2jet);
      RecoRaw_T3jet_angle.push_back(tempJet_T3jet_angle);
      RecoRaw_T3jet      .push_back(tempJet_T3jet);
      RecoRaw_Tpruning   .push_back(tempJet_Tpruning);
      RecoRaw_Ttrimming  .push_back(tempJet_Ttrimming);
      RecoRaw_Taktreclustering.push_back(tempJet_Taktreclustering);
      RecoRaw_Tktreclustering.push_back(tempJet_Tktreclustering);
      RecoRaw_TJet_m1    .push_back(tempJet_TJet_m1);
      RecoRaw_TJet_m2    .push_back(tempJet_TJet_m2);
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
      jetflavor = GetJetTruthFlavor(jettemp, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
      if(debug) cout<<"FillingJet Trimmed: flav="<<jetflavor<<"  pt="<<jettemp.Pt()<<"  m="<<jettemp.M()<<endl;

      if(jetflavor==-1)
        continue;

      /////////////////////////////////
      //Fill variables that will go into ntuple
      /////////////////////////////////
      tempJet_flavor         = jetflavor;
      tempJet_pt             = jettemp.Pt();
      tempJet_eta            = jettemp.Eta();
      tempJet_phi            = jettemp.Phi();
      tempJet_m              = jettemp.M();
      tempJet_Tau21          = GetTau21(groomed_jet);
      tempJet_Tau32          = GetTau32(groomed_jet);
      tempJet_D2             = ecfD2(groomed_jet);

      TSub  T2SubOutputTrim = T_2Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T2jet_angle    = T2SubOutputTrim.min_angle;
      tempJet_T2jet          = T2SubOutputTrim.volatility;

      TSub  T3SubOutputTrim = T_3Subjet(groomed_jet, 0.05, 0.6, 20);
      tempJet_T3jet_angle    = T3SubOutputTrim.min_angle;
      tempJet_T3jet          = T3SubOutputTrim.volatility;

      tempJet_Tpruning       = T_Pruning (groomed_jet, 0.1, 2.0, 20);
      tempJet_Ttrimming      = T_Trimming(groomed_jet, 0.0, 0.1, 20);
      tempJet_Taktreclustering = T_AkTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_Tktreclustering = T_kTreclustering(groomed_jet, 0.05, 0.6, 20);
      tempJet_TJet_m1        = T_Mass(1,groomed_jet);
      tempJet_TJet_m2        = T_Mass(2,groomed_jet);

      RecoRawTrim_flavor     .push_back(tempJet_flavor);
      RecoRawTrim_pt         .push_back(tempJet_pt);
      RecoRawTrim_eta        .push_back(tempJet_eta);
      RecoRawTrim_phi        .push_back(tempJet_phi);
      RecoRawTrim_m          .push_back(tempJet_m);
      RecoRawTrim_Tau21      .push_back(tempJet_Tau21);
      RecoRawTrim_Tau32      .push_back(tempJet_Tau32);
      RecoRawTrim_D2         .push_back(tempJet_D2);
      RecoRawTrim_T2jet_angle.push_back(tempJet_T2jet_angle);
      RecoRawTrim_T2jet      .push_back(tempJet_T2jet);
      RecoRawTrim_T3jet_angle.push_back(tempJet_T3jet_angle);
      RecoRawTrim_T3jet      .push_back(tempJet_T3jet);
      RecoRawTrim_Tpruning   .push_back(tempJet_Tpruning);
      RecoRawTrim_Ttrimming  .push_back(tempJet_Ttrimming);
      RecoRawTrim_Taktreclustering .push_back(tempJet_Taktreclustering);
      RecoRawTrim_Tktreclustering .push_back(tempJet_Tktreclustering);
      RecoRawTrim_TJet_m1    .push_back(tempJet_TJet_m1);
      RecoRawTrim_TJet_m2    .push_back(tempJet_TJet_m2);
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
  TruthRaw_Tau21.clear();
  TruthRaw_Tau32.clear();
  TruthRaw_D2.clear();
  TruthRaw_T2jet_angle.clear();
  TruthRaw_T2jet.clear();
  TruthRaw_T3jet_angle.clear();
  TruthRaw_T3jet.clear();
  TruthRaw_Tpruning.clear();
  TruthRaw_Ttrimming.clear();
  TruthRaw_Taktreclustering.clear();
  TruthRaw_Tktreclustering.clear();
  TruthRaw_TJet_m1.clear();
  TruthRaw_TJet_m2.clear();

  TruthRawTrim_flavor.clear();
  TruthRawTrim_pt.clear();
  TruthRawTrim_eta.clear();
  TruthRawTrim_phi.clear();
  TruthRawTrim_m.clear();
  TruthRawTrim_Tau21.clear();
  TruthRawTrim_Tau32.clear();
  TruthRawTrim_D2.clear();
  TruthRawTrim_T2jet_angle.clear();
  TruthRawTrim_T2jet.clear();
  TruthRawTrim_T3jet_angle.clear();
  TruthRawTrim_T3jet.clear();
  TruthRawTrim_Tpruning.clear();
  TruthRawTrim_Ttrimming.clear();
  TruthRawTrim_Taktreclustering.clear();
  TruthRawTrim_Tktreclustering.clear();
  TruthRawTrim_TJet_m1.clear();
  TruthRawTrim_TJet_m2.clear();

  RecoRaw_flavor.clear();
  RecoRaw_pt.clear();
  RecoRaw_eta.clear();
  RecoRaw_phi.clear();
  RecoRaw_m.clear();
  RecoRaw_Tau21.clear();
  RecoRaw_Tau32.clear();
  RecoRaw_D2.clear();
  RecoRaw_T2jet_angle.clear();
  RecoRaw_T2jet.clear();
  RecoRaw_T3jet_angle.clear();
  RecoRaw_T3jet.clear();
  RecoRaw_Tpruning.clear();
  RecoRaw_Ttrimming.clear();
  RecoRaw_Taktreclustering.clear();
  RecoRaw_Tktreclustering.clear();
  RecoRaw_TJet_m1.clear();
  RecoRaw_TJet_m2.clear();

  RecoRawTrim_flavor.clear();
  RecoRawTrim_pt.clear();
  RecoRawTrim_eta.clear();
  RecoRawTrim_phi.clear();
  RecoRawTrim_m.clear();
  RecoRawTrim_Tau21.clear();
  RecoRawTrim_Tau32.clear();
  RecoRawTrim_D2.clear();
  RecoRawTrim_T2jet_angle.clear();
  RecoRawTrim_T2jet.clear();
  RecoRawTrim_T3jet_angle.clear();
  RecoRawTrim_T3jet.clear();
  RecoRawTrim_Tpruning.clear();
  RecoRawTrim_Ttrimming.clear();
  RecoRawTrim_Taktreclustering.clear();
  RecoRawTrim_Tktreclustering.clear();
  RecoRawTrim_TJet_m1.clear();
  RecoRawTrim_TJet_m2.clear();
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
  for (int p = 0; p < (int)truth_particles.size(); p++) {
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
/// Telescoping Pruning
///=========================================
double T_Pruning(PseudoJet& input, double dcut_min, double dcut_max, int N_dcut) {
  double zcut = 0.1; // single choice of zcut. can be further telescoped
  double Tmass = 0;
  vector<double> ms; ms.clear();
  for (int i = 0; i < N_dcut; i++){
  double dcut = dcut_min + (dcut_max-dcut_min)*i/(N_dcut-1);
  fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, dcut);
  fastjet::PseudoJet pruned_jet = pruner(input);
  Tmass = pruned_jet.m();
  if(Tmass > M0){ms.push_back(Tmass);}
  }
  // getVolatility function provided by TelescopingJets
  if(ms.size()>0) return getVolatility(ms);
  std::cout << "WARNING: zero entries for T_Pruning "<<dcut_min<<" "<<dcut_max<<" "<<N_dcut<<std::endl;
  return -1;
}


///=========================================
/// Telescoping Trimming
///=========================================
double T_Trimming(PseudoJet& input, double fcut_min, double fcut_max, int N_fcut) {
  double Rfilt = 0.1; // single choice of Rfilt. can be further telescoped.
  // used Rfilt = 0.1 for higher pT jets and Rfilt = 0.2 for lower pT jets. 
  double Tmass = 0;
  vector<double> ms; ms.clear();
  for (int i = 0; i < N_fcut; i++){
  double fcut = fcut_min + (fcut_max-fcut_min)*i/(N_fcut-1);
  fastjet::Filter trimmer(Rfilt,fastjet::SelectorPtFractionMin(fcut));
  fastjet::PseudoJet trimmed_jet = trimmer(input);
  Tmass = trimmed_jet.m();
  if(Tmass > M0){ms.push_back(Tmass);}
  }
  // getVolatility function provided by TelescopingJets
  if(ms.size()>0) return getVolatility(ms);
  std::cout << "WARNING: zero entries for T_Trimming "<<fcut_min<<" "<<fcut_max<<" "<<N_fcut<<std::endl;
  return -1;
}


///=========================================
/// Telescoping Anti-kT reclustering
///=========================================
double T_AkTreclustering(PseudoJet& input, double R_min, double R_max, int N_R) {

  vector<double> ms; 

  double Tmass = 0;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
  // used R_min = 0.05 for higher pT jets and R_min = 0.1 for lower pT jets.
  for(int i = 0; i < N_R; i++){

//    R = R_min + i*deltaR;
    fastjet::JetDefinition Tjet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence Tcs(input.constituents(), Tjet_def);

    vector<fastjet::PseudoJet> recoTjets = sorted_by_pt(Tcs.inclusive_jets());
    if(recoTjets.size() < 1) {
        std::cout <<"Warning: recluster number of subjet is "<< recoTjets.size() << std::endl;
        continue;
    }
    if(recoTjets.size() == 1){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tmass = Tsubjet1.M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    else if(recoTjets.size() >= 2){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tsubjet2.SetPxPyPzE(recoTjets[1].px(),recoTjets[1].py(),recoTjets[1].pz(),recoTjets[1].e());
        Tmass = (Tsubjet1+Tsubjet2).M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    R += deltaR;
  }
  // getVolatility function provided by TelescopingJets
  if(ms.size()>0) return getVolatility(ms);
  std::cout << "WARNING: zero entries for T_AkTreclustering "<<R_min<<" "<<R_max<<" "<<N_R<<std::endl;
  return -1;
}


///=========================================
/// Telescoping kT reclustering
///=========================================
double T_kTreclustering(PseudoJet& input, double R_min, double R_max, int N_R) {

  vector<double> ms; 

  double Tmass = 0;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
  for(int i = 0; i < N_R; i++){

//    R = R_min + i*deltaR;
    fastjet::JetDefinition Tjet_def(fastjet::kt_algorithm, R);
    fastjet::ClusterSequence Tcs(input.constituents(), Tjet_def);

    vector<fastjet::PseudoJet> recoTjets = sorted_by_pt(Tcs.inclusive_jets());
    if(recoTjets.size() < 1) {
        std::cout <<"Warning: recluster number of subjet is "<< recoTjets.size() << std::endl;
        continue;
    }
    if(recoTjets.size() == 1){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tmass = Tsubjet1.M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    else if(recoTjets.size() >= 2){
        Tsubjet1.SetPxPyPzE(recoTjets[0].px(),recoTjets[0].py(),recoTjets[0].pz(),recoTjets[0].e());
        Tsubjet2.SetPxPyPzE(recoTjets[1].px(),recoTjets[1].py(),recoTjets[1].pz(),recoTjets[1].e());
        Tmass = (Tsubjet1+Tsubjet2).M();
	if(Tmass > M0){ms.push_back(Tmass);}
    }
    R += deltaR;     
  }
  // getVolatility function provided by TelescopingJets
  if(ms.size()>0) return getVolatility(ms);
  std::cout << "WARNING: zero entries for T_kTreclustering "<<R_min<<" "<<R_max<<" "<<N_R<<std::endl;
  return -1;
}


///=========================================
/// Telescoping Subjet
///=========================================

TSub T_2Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
  vector<double> m12s;
  double beta = 1.0;
  double Tmass = 0;
  fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
  fastjet::contrib::Nsubjettiness nSub(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
  double tau2 = nSub.result(input);
  std::vector<fastjet::PseudoJet> tau2axes = nSub.currentAxes();
  tau_axis1.SetPxPyPzE(tau2axes[0].px(),tau2axes[0].py(),tau2axes[0].pz(),tau2axes[0].e());
  tau_axis2.SetPxPyPzE(tau2axes[1].px(),tau2axes[1].py(),tau2axes[1].pz(),tau2axes[1].e());
  double D_min = tau_axis1.DeltaR(tau_axis2);

  double d1,d2;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;
 
  for (int i = 0; i < N_R; i++){
//    R = R_min + i*deltaR;
    Tsubjet1.SetPxPyPzE(0,0,0,0);
    Tsubjet2.SetPxPyPzE(0,0,0,0);

    for (unsigned int c = 0; c < input.constituents().size(); c++) {
      particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
      d1 = particle.DeltaR(tau_axis1);
      d2 = particle.DeltaR(tau_axis2);
      if (d1 <= R && d1 < d2 ){
        Tsubjet1 = Tsubjet1 + particle;
      }
      else if (d2 <= R && d2 < d1 ){
        Tsubjet2 = Tsubjet2 + particle;
      }
    }
    Tmass = (Tsubjet1 + Tsubjet2).M();
    if(Tmass > M0){m12s.push_back(Tmass);}
    R += deltaR;
  }
  TSub result;
  result.min_angle = D_min;
  if(m12s.size()>0) result.volatility = getVolatility(m12s);
  else {
    result.volatility = -1;
    std::cout << "WARNING: zero entries for T_2Subjet "<<R_min<<" "<<R_max<<" "<<N_R<<std::endl;
  }
  return result;
}



TSub T_3Subjet(PseudoJet& input, double R_min, double R_max, int N_R){
  vector<double> m123s; m123s.clear();
  // m123 is the invariant mass of the three subjets
  double beta = 1.0;
  double Tmass = 0;
  fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
  fastjet::contrib::Nsubjettiness nSub(3, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
  double tau3 = nSub.result(input);
  std::vector<fastjet::PseudoJet> tau3axes = nSub.currentAxes();
  tau_axis1.SetPxPyPzE(tau3axes[0].px(),tau3axes[0].py(),tau3axes[0].pz(),tau3axes[0].e());
  tau_axis2.SetPxPyPzE(tau3axes[1].px(),tau3axes[1].py(),tau3axes[1].pz(),tau3axes[1].e());
  tau_axis3.SetPxPyPzE(tau3axes[2].px(),tau3axes[2].py(),tau3axes[2].pz(),tau3axes[2].e());
  double tau3_D12 = tau_axis1.DeltaR(tau_axis2);
  double tau3_D13 = tau_axis1.DeltaR(tau_axis3);
  double tau3_D23 = tau_axis2.DeltaR(tau_axis3);
  double D_min = tau3_D12;
  if(tau3_D13 < D_min ) D_min = tau3_D13;
  if(tau3_D23 < D_min ) D_min = tau3_D23;

  double d1, d2, d3;
  double deltaR = (R_max - R_min)/(N_R-1);
  double R      = R_min;

  for (int i = 0; i < N_R; i++){
//    R = R_min + i*deltaR;
    Tsubjet1.SetPxPyPzE(0,0,0,0);
    Tsubjet2.SetPxPyPzE(0,0,0,0);
    Tsubjet3.SetPxPyPzE(0,0,0,0);

    for (unsigned int c = 0; c < input.constituents().size(); c++) {
      particle.SetPxPyPzE(input.constituents()[c].px(),input.constituents()[c].py(),input.constituents()[c].pz(),input.constituents()[c].e());
      d1 = particle.DeltaR(tau_axis1);
      d2 = particle.DeltaR(tau_axis2);
      d3 = particle.DeltaR(tau_axis3);
      if (d1 <= R && d1 < d2 && d1 < d3){
        Tsubjet1 = Tsubjet1 + particle;
      }
      else if (d2 <= R && d2 < d1 && d2 < d3){
        Tsubjet2 = Tsubjet2 + particle;
      }
      else if (d3 <= R && d3 < d1 && d3 < d2){
        Tsubjet3 = Tsubjet3 + particle;
      }
    }
    Tmass = (Tsubjet1 + Tsubjet2 + Tsubjet3).M();
    if(Tmass > M0){m123s.push_back(Tmass);}
    R += deltaR;
  }
  TSub result;
  result.min_angle = D_min;
  if(m123s.size()>0) result.volatility = getVolatility(m123s);
  else {
    result.volatility = -1;
    std::cout << "WARNING: zero entries for T_3Subjet "<<R_min<<" "<<R_max<<" "<<N_R<<std::endl;
  }
  return result;
}


///=========================================
/// Telescoping N-subjettiness
///=========================================
double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> taus; taus.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(beta);
    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
//    fastjet::contrib::Nsubjettiness nsub(N, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
    taus.push_back(nsub(input));
  }
  // getVolatility function provided by TelescopingJets
  if(taus.size()>0) return getVolatility(taus);
  std::cout << "WARNING: zero entries for T_Nsubjettiness "<<beta_min<<" "<<beta_max<<" "<<N_beta<<std::endl;
  return -1;
}

double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> taus; taus.clear();
  for (int i = 0; i < N_beta; i++) {

    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);

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
  if(taus.size()>0) return getVolatility(taus);
  std::cout << "WARNING: zero entries for T_NsubjettinessRatio "<<beta_min<<" "<<beta_max<<" "<<N_beta<<std::endl;
  return -1;
}


///=========================================
/// Telescoping Energy Correlators
///=========================================
double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorC2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  if(ecfs.size()>0) return getVolatility(ecfs);
  std::cout << "WARNING: zero entries for T_EnergyCorrelator_C2 "<<beta_min<<" "<<beta_max<<" "<<N_beta<<std::endl;
  return -1;
}

double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorD2 ecf(beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  if(ecfs.size()>0) return getVolatility(ecfs);
  std::cout << "WARNING: zero entries for T_EnergyCorrelator_D2 "<<beta_min<<" "<<beta_max<<" "<<N_beta<<std::endl;
  return -1;
}

double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max, int N_beta) {
  vector<double> ecfs; ecfs.clear();
  for (int i = 0; i < N_beta; i++) {
    double beta = beta_min + i*(beta_max - beta_min)/(N_beta-1);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecf(3, beta);
    ecfs.push_back(ecf(input));
  }
  // getVolatility function provided by TelescopingJets
  if(ecfs.size()>0) return getVolatility(ecfs);
  std::cout << "WARNING: zero entries for T_EnergyCorrelator_C3 "<<beta_min<<" "<<beta_max<<" "<<N_beta<<std::endl;
  return -1;
}

///========================================
int GetJetTruthFlavor(TLorentzVector jettemp,
                      TLorentzVector truth_t1,
                      TLorentzVector truth_t2,
                      TLorentzVector truth_W1,
                      TLorentzVector truth_W2,
                      TLorentzVector truth_H,
                      int debug){
  if(debug){
    cout<<"DeltaR:   "<<endl
        <<"dRMatch:  "<<dR_match<<endl
        <<"q1:       "<<jettemp.DeltaR(truth_q1)<<endl
        <<"q2:       "<<jettemp.DeltaR(truth_q2)<<endl
        <<"W1:       "<<jettemp.DeltaR(truth_W1)<<endl
        <<"W2:       "<<jettemp.DeltaR(truth_W2)<<endl
        <<"H:        "<<jettemp.DeltaR(truth_H)<<endl
        <<"t1:       "<<jettemp.DeltaR(truth_t1)<<endl
        <<"t2:       "<<jettemp.DeltaR(truth_t2)<<endl;
  }
  int jetflavor = -1;
  if     (jettemp.DeltaR(truth_q1)<dR_match || jettemp.DeltaR(truth_q2)<dR_match){
    jetflavor = 0;
  }
  else if(jettemp.DeltaR(truth_W1)<dR_match || jettemp.DeltaR(truth_W2)<dR_match){
    jetflavor = 1;
  }
  else if(jettemp.DeltaR(truth_t1)<dR_match || jettemp.DeltaR(truth_t2)<dR_match){
    jetflavor = 3;
  }
  else if(jettemp.DeltaR(truth_H)<dR_match){
    jetflavor = 3;
  }
  else{
    jetflavor = -1;
  }

  if(debug) cout<<"Found jet truth flavor: "<<jetflavor<<endl;

  return jetflavor;
}


double GetTau21(PseudoJet& input){

  float tau21=-1;

  //N-subjettiness
  fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
  fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
//  fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
//  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);

  float tau1 = nsub1(input);
  float tau2 = nsub2(input);

  if(tau1>0)
    tau21 = tau2/tau1;

  return tau21;

}

double GetTau32(PseudoJet& input){

  float tau32=-1;

  //N-subjettiness
  fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
  fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::OnePass_KT_Axes(), nsubMeasure);
//  fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);
//  fastjet::contrib::Nsubjettiness nsub3(3, fastjet::contrib::WTA_KT_Axes(), nsubMeasure);

  float tau2 = nsub2(input);
  float tau3 = nsub3(input);

  if(tau2>0)
    tau32 = tau3/tau2;

  return tau32;

}


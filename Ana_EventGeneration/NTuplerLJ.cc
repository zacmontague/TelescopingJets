/*----------------------------------------------------------------------

TITLE : NTuplerLJ.cc

DESCRIPTION : Takes input final state particle ntuple from Ana_EventGeneration
and outputs a flat ntuple that contains jet properties created via fastjet
in this code to be used for post-analysis using Ana_MiniNTupleAnalysis. (NOTE!!)
All of the TelescopingJets code from the fastjet contrib should be contained within this code

COMPILE :
$ source compile.sh

RUN :
$ ./NTuplerLJ <type> <input> <output>

type   : 0 = dijet , 1 = G*->W+W- , 2 = ttbar
input  : Input file from Ana_MiniNTupleAnalysis
output : Anything you want - but being logical

//----------------------------------------------------------------------*/


#include "NTuplerLJ.h"

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

  treeout->Branch("TruthRaw_flavor",        &truth_moments.flavor);
  treeout->Branch("TruthRaw_pt",            &truth_moments.pt);
  treeout->Branch("TruthRaw_eta",           &truth_moments.eta);
  treeout->Branch("TruthRaw_phi",           &truth_moments.phi);
  treeout->Branch("TruthRaw_m",             &truth_moments.m);
  treeout->Branch("TruthRaw_Tau21",         &truth_moments.Tau21);
  treeout->Branch("TruthRaw_Tau32",         &truth_moments.Tau32);
  treeout->Branch("TruthRaw_D2",            &truth_moments.D2);
  treeout->Branch("TruthRaw_T2jet_angle",   &truth_moments.T2jet_angle);
  treeout->Branch("TruthRaw_T2jet",         &truth_moments.T2jet);
  treeout->Branch("TruthRaw_T3jet_angle",   &truth_moments.T3jet_angle);
  treeout->Branch("TruthRaw_T3jet",         &truth_moments.T3jet);
  treeout->Branch("TruthRaw_Tpruning",      &truth_moments.Tpruning);
  treeout->Branch("TruthRaw_Ttrimming",     &truth_moments.Ttrimming);
  treeout->Branch("TruthRaw_Taktreclustering",	&truth_moments.Taktreclustering);
  treeout->Branch("TruthRaw_Tktreclustering",	&truth_moments.Tktreclustering);
  treeout->Branch("TruthRaw_TJet_m1",       &truth_moments.TJet_m1);
  treeout->Branch("TruthRaw_TJet_m2",       &truth_moments.TJet_m2);

  treeout->Branch("TruthRawTrimflavor",        &truthtrim_moments.flavor);
  treeout->Branch("TruthRawTrimpt",            &truthtrim_moments.pt);
  treeout->Branch("TruthRawTrimeta",           &truthtrim_moments.eta);
  treeout->Branch("TruthRawTrimphi",           &truthtrim_moments.phi);
  treeout->Branch("TruthRawTrimm",             &truthtrim_moments.m);
  treeout->Branch("TruthRawTrimTau21",         &truthtrim_moments.Tau21);
  treeout->Branch("TruthRawTrimTau32",         &truthtrim_moments.Tau32);
  treeout->Branch("TruthRawTrimD2",            &truthtrim_moments.D2);
  treeout->Branch("TruthRawTrimT2jet_angle",   &truthtrim_moments.T2jet_angle);
  treeout->Branch("TruthRawTrimT2jet",         &truthtrim_moments.T2jet);
  treeout->Branch("TruthRawTrimT3jet_angle",   &truthtrim_moments.T3jet_angle);
  treeout->Branch("TruthRawTrimT3jet",         &truthtrim_moments.T3jet);
  treeout->Branch("TruthRawTrimTpruning",      &truthtrim_moments.Tpruning);
  treeout->Branch("TruthRawTrimTtrimming",     &truthtrim_moments.Ttrimming);
  treeout->Branch("TruthRawTrimTaktreclustering",	&truthtrim_moments.Taktreclustering);
  treeout->Branch("TruthRawTrimTktreclustering",	&truthtrim_moments.Tktreclustering);
  treeout->Branch("TruthRawTrimTJet_m1",       &truthtrim_moments.TJet_m1);
  treeout->Branch("TruthRawTrimTJet_m2",       &truthtrim_moments.TJet_m2);

  treeout->Branch("RecoRawflavor",        &reco_moments.flavor);
  treeout->Branch("RecoRawpt",            &reco_moments.pt);
  treeout->Branch("RecoRaweta",           &reco_moments.eta);
  treeout->Branch("RecoRawphi",           &reco_moments.phi);
  treeout->Branch("RecoRawm",             &reco_moments.m);
  treeout->Branch("RecoRawTau21",         &reco_moments.Tau21);
  treeout->Branch("RecoRawTau32",         &reco_moments.Tau32);
  treeout->Branch("RecoRawD2",            &reco_moments.D2);
  treeout->Branch("RecoRawT2jet_angle",   &reco_moments.T2jet_angle);
  treeout->Branch("RecoRawT2jet",         &reco_moments.T2jet);
  treeout->Branch("RecoRawT3jet_angle",   &reco_moments.T3jet_angle);
  treeout->Branch("RecoRawT3jet",         &reco_moments.T3jet);
  treeout->Branch("RecoRawTpruning",      &reco_moments.Tpruning);
  treeout->Branch("RecoRawTtrimming",     &reco_moments.Ttrimming);
  treeout->Branch("RecoRawTaktreclustering",	&reco_moments.Taktreclustering);
  treeout->Branch("RecoRawTktreclustering",	&reco_moments.Tktreclustering);
  treeout->Branch("RecoRawTJet_m1",       &reco_moments.TJet_m1);
  treeout->Branch("RecoRawTJet_m2",       &reco_moments.TJet_m2);

  treeout->Branch("RecoRawTrimflavor",        &recotrim_moments.flavor);
  treeout->Branch("RecoRawTrimpt",            &recotrim_moments.pt);
  treeout->Branch("RecoRawTrimeta",           &recotrim_moments.eta);
  treeout->Branch("RecoRawTrimphi",           &recotrim_moments.phi);
  treeout->Branch("RecoRawTrimm",             &recotrim_moments.m);
  treeout->Branch("RecoRawTrimTau21",         &recotrim_moments.Tau21);
  treeout->Branch("RecoRawTrimTau32",         &recotrim_moments.Tau32);
  treeout->Branch("RecoRawTrimD2",            &recotrim_moments.D2);
  treeout->Branch("RecoRawTrimT2jet_angle",   &recotrim_moments.T2jet_angle);
  treeout->Branch("RecoRawTrimT2jet",         &recotrim_moments.T2jet);
  treeout->Branch("RecoRawTrimT3jet_angle",   &recotrim_moments.T3jet_angle);
  treeout->Branch("RecoRawTrimT3jet",         &recotrim_moments.T3jet);
  treeout->Branch("RecoRawTrimTpruning",      &recotrim_moments.Tpruning);
  treeout->Branch("RecoRawTrimTtrimming",     &recotrim_moments.Ttrimming);
  treeout->Branch("RecoRawTrimTaktreclustering",	&recotrim_moments.Taktreclustering);
  treeout->Branch("RecoRawTrimTktreclustering",	&recotrim_moments.Tktreclustering);
  treeout->Branch("RecoRawTrimTJet_m1",       &recotrim_moments.TJet_m1);
  treeout->Branch("RecoRawTrimTJet_m2",       &recotrim_moments.TJet_m2);

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
  //for (Long64_t jentry=0; jentry<4000/*nEvents*/; jentry++) {
  //for (Long64_t jentry=4000; jentry<8000/*nEvents*/; jentry++) {
  //for (Long64_t jentry=8000; jentry<12000/*nEvents*/; jentry++) {
  //for (Long64_t jentry=12000; jentry<16000/*nEvents*/; jentry++) {
  //for (Long64_t jentry=16000; jentry<nEvents; jentry++) {
  //for (Long64_t jentry=12000; jentry<15000/*nEvents*/; jentry++) {

    if(jentry%100==0)
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


    /////////////////////////////////////////////
    //make pseudocalorimeter cells
    //////////////////////////////////////////////
    vector<PseudoJet> calo_cells        = ToyCalorimeter(input_particles);
    //vector<PseudoJet> calo_cells_Pileup = ToyCalorimeter(input_particles_Pileup);

  
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
  
    //////////////////////////////////////////////
    // get the jets ordered in pt
    //////////////////////////////////////////////
    fastjet::ClusterSequence clust_seq_truth(input_particles, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_truth = sorted_by_pt(clust_seq_truth.inclusive_jets(5.0));
    vector<fastjet::PseudoJet> trimmed_jets_truth = vector<fastjet::PseudoJet>();
    for (vector<PseudoJet>::iterator jiter_truth=inclusive_jets_truth.begin();
          jiter_truth!=inclusive_jets_truth.end(); jiter_truth++) {
      trimmed_jets_truth.push_back(f(*jiter_truth));
    }
    trimmed_jets_truth = sorted_by_pt(trimmed_jets_truth);

    fastjet::ClusterSequence clust_seq_reco(input_particles, jet_def);
    vector<fastjet::PseudoJet> inclusive_jets_reco = sorted_by_pt(clust_seq_reco.inclusive_jets(5.0));
    vector<fastjet::PseudoJet> trimmed_jets_reco = vector<fastjet::PseudoJet>();
    for (vector<PseudoJet>::iterator jiter_reco=inclusive_jets_reco.begin();
          jiter_reco!=inclusive_jets_reco.end(); jiter_reco++) {
      trimmed_jets_reco.push_back(f(*jiter_reco));
    }
    trimmed_jets_reco = sorted_by_pt(trimmed_jets_reco);

    if(debug){
      // label the columns
      cout<<"jet#  pt  eta  phi  mass"<<endl;
      cout<<"Inclusive"<<endl;
      // print out the details for each jet
      for (unsigned int i = 0; i < inclusive_jets_truth.size(); i++) {
        cout<<i<<"  "<<inclusive_jets_truth[i].pt()
               <<"  "<<inclusive_jets_truth[i].eta()
               <<"  "<<inclusive_jets_truth[i].phi()
               <<"  "<<inclusive_jets_truth[i].m()<<endl;
      }
    }
       
    /////////////////////////////////////////////
    // fill tree
    //////////////////////////////////////////////
    truth_moments.SetAll(GetMoments(inclusive_jets_truth, NULL));
    truthtrim_moments.SetAll(GetMoments(trimmed_jets_truth, NULL));
    reco_moments.SetAll(GetMoments(inclusive_jets_reco, NULL));
    recotrim_moments.SetAll(GetMoments(trimmed_jets_reco, NULL));

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

//////////////////////////////////////////////
// Take the leading jet matched to truth_ljet
// and calculate all the variables for the ntuple
//////////////////////////////////////////////
JetMoments *GetMoments(vector<PseudoJet> inclusive_jets, PseudoJet *truth_ljet) {
  bool debug = false;
  
  /////////////////////////////
  // Take leading jet matched to truth_ljet
  /////////////////////////////
  if(debug) cout<<"Processing jet"<<endl;
  PseudoJet jettemp;
  bool matched = false;
  if (truth_ljet == NULL) {
    jettemp = inclusive_jets.at(0);
    matched = true;
  } else {
    for(int ijet=0; ijet<inclusive_jets.size(); ijet++){
      if (truth_ljet->delta_R(inclusive_jets.at(ijet)) < 1.0) {
        matched = true;
        break;
      }
    }
  }
  if (!matched) {
    cout<<"No matching jet found. Skipping "<<endl;
    return NULL;
  }

  return CalcMoments(jettemp);
}

JetMoments *CalcMoments(PseudoJet &jet) {
  bool debug = false;
  PseudoJet jettemp = jet;

  /////////////////////////////////
  //Getting truth label for filling into ntuple
  /////////////////////////////////
  TLorentzVector jettempLV;
  jettempLV.SetPtEtaPhiM(jettemp.pt(),
                         jettemp.eta(),
                         jettemp.phi(),
                         jettemp.m());
  jetflavor = GetJetTruthFlavor(jettempLV, truth_t1, truth_t2, truth_W1, truth_W2, truth_H, debug);
  if(debug) cout<<"FillingJet Raw   : flav="<<jetflavor<<"  pt="<<jettemp.pt()<<"  m="<<jettemp.m()<<endl;

  if(jetflavor==-1) {
    cout<<"Matching jet is wrong flavor. Skipping"<<endl; 
    return NULL;
  }

  //////////////////////////////////////////////
  //Set up Telescoping jets (this looks like the Telescoping reclustering)
  //////////////////////////////////////////////
  fastjet::contrib::KT_Axes axes_def;
  std::vector<double> r_values;
  int N_r = 20;
  double r_min = 0.1;
  double r_max = 0.6;
  for(int i=0; i < N_r; i++){
    r_values.push_back( r_min+i*(r_max-r_min)/(N_r-1) );
  }
  TelescopingJets T_Mass(axes_def,r_values);

  /////////////////////////////////
  //Calculate variables that will go into ntuple
  /////////////////////////////////
  if (debug) cout<<"Ready for calculation"<<endl; 
  tempJet_flavor         = jetflavor;
  tempJet_pt             = jettemp.pt();
  tempJet_eta            = jettemp.eta();
  tempJet_phi            = jettemp.phi();
  tempJet_m              = jettemp.m();
  tempJet_Tau21          = GetTau21(jettemp);
  tempJet_Tau32          = GetTau32(jettemp);
  tempJet_D2             = ecfD2(jettemp);

  TSub  T2SubOutput     = T_2Subjet(jettemp, 0.05, 0.6, 20);
  tempJet_T2jet_angle    = T2SubOutput.min_angle;
  tempJet_T2jet          = T2SubOutput.volatility;

  TSub  T3SubOutput     = T_3Subjet(jettemp, 0.05, 0.6, 20);
  tempJet_T3jet_angle    = T3SubOutput.min_angle;
  tempJet_T3jet          = T3SubOutput.volatility;

  tempJet_Tpruning       = T_Pruning (jettemp, 0.1, 2.0, 20);
  tempJet_Ttrimming      = T_Trimming(jettemp, 0.0, 0.1, 20);
  tempJet_Taktreclustering = T_AkTreclustering(jettemp, 0.05, 0.6, 20);
  tempJet_Tktreclustering = T_kTreclustering(jettemp, 0.05, 0.6, 20);
  tempJet_TJet_m1        = T_Mass(1,jettemp);
  tempJet_TJet_m2        = T_Mass(2,jettemp);

  /////////////////////////////////
  //Fill variables for output
  /////////////////////////////////
  if (debug) cout<<"Storing values"<<endl; 
  JetMoments *output = new JetMoments();
  output->flavor = tempJet_flavor;
  output->pt = tempJet_pt;
  output->eta = tempJet_eta;
  output->phi = tempJet_phi;
  output->m = tempJet_m;
  output->Tau21 = tempJet_Tau21;
  output->Tau32 = tempJet_Tau32;
  output->D2 = tempJet_D2;
  output->T2jet_angle = tempJet_T2jet_angle;
  output->T2jet = tempJet_T2jet;
  output->T3jet_angle = tempJet_T3jet_angle;
  output->T3jet = tempJet_T3jet;
  output->Tpruning = tempJet_Tpruning;
  output->Ttrimming = tempJet_Ttrimming;
  output->Taktreclustering = tempJet_Taktreclustering;
  output->Tktreclustering = tempJet_Tktreclustering;
  output->TJet_m1 = tempJet_TJet_m1;
  output->TJet_m2 = tempJet_TJet_m2;

  return output;
}


///=========================================
/// Reset Branches
///=========================================
void ResetBranches(){

  NumberOfVertices = 0;

  truth_moments.flavor = 0;
  truth_moments.pt = 0;
  truth_moments.eta = 0;
  truth_moments.phi = 0;
  truth_moments.m = 0;
  truth_moments.Tau21 = 0;
  truth_moments.Tau32 = 0;
  truth_moments.D2 = 0;
  truth_moments.T2jet_angle = 0;
  truth_moments.T2jet = 0;
  truth_moments.T3jet_angle = 0;
  truth_moments.T3jet = 0;
  truth_moments.Tpruning = 0;
  truth_moments.Ttrimming = 0;
  truth_moments.Taktreclustering = 0;
  truth_moments.Tktreclustering = 0;
  truth_moments.TJet_m1 = 0;
  truth_moments.TJet_m2 = 0;

  truthtrim_moments.flavor = 0;
  truthtrim_moments.pt = 0;
  truthtrim_moments.eta = 0;
  truthtrim_moments.phi = 0;
  truthtrim_moments.m = 0;
  truthtrim_moments.Tau21 = 0;
  truthtrim_moments.Tau32 = 0;
  truthtrim_moments.D2 = 0;
  truthtrim_moments.T2jet_angle = 0;
  truthtrim_moments.T2jet = 0;
  truthtrim_moments.T3jet_angle = 0;
  truthtrim_moments.T3jet = 0;
  truthtrim_moments.Tpruning = 0;
  truthtrim_moments.Ttrimming = 0;
  truthtrim_moments.Taktreclustering = 0;
  truthtrim_moments.Tktreclustering = 0;
  truthtrim_moments.TJet_m1 = 0;
  truthtrim_moments.TJet_m2 = 0;

  reco_moments.flavor = 0;
  reco_moments.pt = 0;
  reco_moments.eta = 0;
  reco_moments.phi = 0;
  reco_moments.m = 0;
  reco_moments.Tau21 = 0;
  reco_moments.Tau32 = 0;
  reco_moments.D2 = 0;
  reco_moments.T2jet_angle = 0;
  reco_moments.T2jet = 0;
  reco_moments.T3jet_angle = 0;
  reco_moments.T3jet = 0;
  reco_moments.Tpruning = 0;
  reco_moments.Ttrimming = 0;
  reco_moments.Taktreclustering = 0;
  reco_moments.Tktreclustering = 0;
  reco_moments.TJet_m1 = 0;
  reco_moments.TJet_m2 = 0;

  recotrim_moments.flavor = 0;
  recotrim_moments.pt = 0;
  recotrim_moments.eta = 0;
  recotrim_moments.phi = 0;
  recotrim_moments.m = 0;
  recotrim_moments.Tau21 = 0;
  recotrim_moments.Tau32 = 0;
  recotrim_moments.D2 = 0;
  recotrim_moments.T2jet_angle = 0;
  recotrim_moments.T2jet = 0;
  recotrim_moments.T3jet_angle = 0;
  recotrim_moments.T3jet = 0;
  recotrim_moments.Tpruning = 0;
  recotrim_moments.Ttrimming = 0;
  recotrim_moments.Taktreclustering = 0;
  recotrim_moments.Tktreclustering = 0;
  recotrim_moments.TJet_m1 = 0;
  recotrim_moments.TJet_m2 = 0;
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


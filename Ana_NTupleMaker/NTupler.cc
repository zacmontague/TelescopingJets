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

int main(){

  //debug flag
  bool debug=false;

  //////////////////////////////////////////////
  //INPUT
  //////////////////////////////////////////////
  //get input file and tree
  TString infilename="../Ana_EventGeneration/output_WprimeWZ.root";
  filein = new TFile( infilename );
  treein = (TTree*)filein->Get( "tree" );
  treein->Print();

  //set up branch linking to addresses
  treein->SetBranchAddress("fspart_id", &fspart_id);
  treein->SetBranchAddress("fspart_pt", &fspart_pt);
  treein->SetBranchAddress("fspart_eta",&fspart_eta);
  treein->SetBranchAddress("fspart_phi",&fspart_phi);
  treein->SetBranchAddress("fspart_m",  &fspart_m);

  treein->SetBranchAddress("truth_quark1_pt",  &truth_quark1_pt);
  treein->SetBranchAddress("truth_quark1_eta", &truth_quark1_eta);
  treein->SetBranchAddress("truth_quark1_phi", &truth_quark1_phi);
  treein->SetBranchAddress("truth_quark1_m",   &truth_quark1_m);

  treein->SetBranchAddress("truth_quark2_pt",  &truth_quark2_pt);
  treein->SetBranchAddress("truth_quark2_eta", &truth_quark2_eta);
  treein->SetBranchAddress("truth_quark2_phi", &truth_quark2_phi);
  treein->SetBranchAddress("truth_quark2_m",   &truth_quark2_m);

  treein->SetBranchAddress("truth_top1_pt",  &truth_top1_pt);
  treein->SetBranchAddress("truth_top1_eta", &truth_top1_eta);
  treein->SetBranchAddress("truth_top1_phi", &truth_top1_phi);
  treein->SetBranchAddress("truth_top1_m",   &truth_top1_m);

  treein->SetBranchAddress("truth_top2_pt",  &truth_top2_pt);
  treein->SetBranchAddress("truth_top2_eta", &truth_top2_eta);
  treein->SetBranchAddress("truth_top2_phi", &truth_top2_phi);
  treein->SetBranchAddress("truth_top2_m",   &truth_top2_m);

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
  treeout = new TTree("treeout","treeout");

  treeout->Branch("Truth_flavor", &Truth_flavor);
  treeout->Branch("Truth_pt",     &Truth_pt);
  treeout->Branch("Truth_eta",    &Truth_eta);
  treeout->Branch("Truth_phi",    &Truth_phi);
  treeout->Branch("Truth_m",      &Truth_m);
  treeout->Branch("Truth_tjet1",  &Truth_tjet1);
  treeout->Branch("Truth_tjet2",  &Truth_tjet2);


//   treeout->Branch("TruthAK10Trim_pt",  &TruthAK10Trim_pt);
//   treeout->Branch("TruthAK10Trim_eta", &TruthAK10Trim_eta);
//   treeout->Branch("TruthAK10Trim_phi", &TruthAK10Trim_phi);
//   treeout->Branch("TruthAK10Trim_m",   &TruthAK10Trim_m);

  //////////////////////////////////////////////
  //main event loop
  //////////////////////////////////////////////
  nEvents = treein->GetEntries();
  cout<<"Number of events: "<<nEvents<<endl;

  for (Long64_t jentry=0; jentry<nEvents; jentry++) {

    if(debug)
      if(jentry>1)
        continue;

    cout<<jentry<<endl;

    filein->cd();
    treein->GetEntry(jentry);

    if(debug) cout<<"TruthTop: "<<truth_top1_m<<"  "<<truth_top2_m<<endl;

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
    //TODO

    //////////////////////////////////////////////
    //make pseudocalorimeter cells
    //////////////////////////////////////////////
    //TODO
    vector<PseudoJet> calo_cells = ToyCalorimeter(input_particles);

    //////////////////////////////////////////////
    // get the resulting jets ordered in pt
    //////////////////////////////////////////////
    double R = 1.0;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    double ptmin = 5.0;
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

    // label the columns
    cout<<"jet#  pt  eta  phi  mass"<<endl;
    cout<<"Inclusive"<<endl;
    // print out the details for each jet
    for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
      cout<<i<<"  "<<inclusive_jets[i].pt()<<"  "<<inclusive_jets[i].eta()<<"  "<<inclusive_jets[i].phi()<<"  "<<inclusive_jets[i].m()<<endl;
    }



    //////////////////////////////////////////////
    //Testing the Telescoping Jets contrib - setup done here
    //////////////////////////////////////////////
    fastjet::contrib::KT_Axes axes_def;
    std::vector<double> r_values;
    for(int i=0; i < 20; i++){
      r_values.push_back( 0.1+i*(0.6-0.1)/(20-1) );
    }
    TelescopingJets tjet(axes_def,r_values);

    //////////////////////////////////////////////
    //Setup calculations for other algorithms
    //////////////////////////////////////////////
    //N-subjettiness
    fastjet::contrib::UnnormalizedMeasure nsubMeasure(1.);
    fastjet::contrib::Nsubjettiness nsub1(1, fastjet::contrib::WTA_KT_Axes, nsubMeasure);
    fastjet::contrib::Nsubjettiness nsub2(2, fastjet::contrib::WTA_KT_Axes, nsubMeasure);

    //Energy correlation functions
    fastjet::contrib::EnergyCorrelatorDoubleRatio ecfC2(2, 1.);
    fastjet::contrib::EnergyCorrelatorD2 ecfD2(1.);

    //////////////////////////////////////////////
    // Filtering with a pt cut as for trimming (arXiv:0912.1342)
    //////////////////////////////////////////////
//     vector<PseudoJet> inclusive_jets_trimmed;
//     inclusive_jets_trimmed.clear();
//     double Rtrim = 0.2;
//     double ptfrac = 0.05;
//     Transformer *trimmer = new Filter(JetDefinition(kt_algorithm, Rtrim), SelectorPtFractionMin(ptfrac) );
//     const Transformer & f = *trimmer;
//     //groom jets
//     printf("Trimmed");
//     for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
//       PseudoJet groomed_jet = f(inclusive_jets[i]);
//       inclusive_jets_trimmed.push_back(groomed_jet);
//       cout<<i<<"  "<<inclusive_jets_trimmed[i].pt()<<"  "<<inclusive_jets_trimmed[i].eta()<<"  "<<inclusive_jets_trimmed[i].phi()<<"  "<<inclusive_jets_trimmed[i].m()<<endl;
//     }


    /////////////////////////////////////////////
    //Get truth objects for truth matching
    /////////////////////////////////////////////
    TLorentzVector truth_W;
    truth_W.SetPtEtaPhiM(truth_W_pt,truth_W_eta,truth_W_phi,truth_W_m);
    TLorentzVector truth_Z;
    truth_Z.SetPtEtaPhiM(truth_Z_pt,truth_Z_eta,truth_Z_phi,truth_Z_m);


    //fill output ntuple
    Truth_pt  = 64.0;
    Truth_eta = 64.0;
    Truth_phi = 64.0;
    Truth_m   = 64.0;

    if(inclusive_jets.size()>0){
      Truth_pt  = inclusive_jets.at(0).pt();
      Truth_eta = inclusive_jets.at(0).eta();
      Truth_phi = inclusive_jets.at(0).phi();
      Truth_m   = inclusive_jets.at(0).m();
      treeout->Fill();
    }


    for(int ijet=0; ijet<inclusive_jets.size(); ijet++){
      TLorentzVector jettemp;
      jettemp.SetPtEtaPhiM(inclusive_jets.at(ijet).pt(),inclusive_jets.at(ijet).eta(),inclusive_jets.at(ijet).phi(),inclusive_jets.at(ijet).m());

      /////////////////////////////////
      //Getting truth label for filling into ntuple
      /////////////////////////////////
      cout<<"DeltaR: "<<endl
          <<"W:  "<<jettemp.DeltaR(truth_W)<<endl
          <<"Z:  "<<jettemp.DeltaR(truth_Z)<<endl;

      jetflavor = -1;
      if(jettemp.DeltaR(truth_W)<0.8){
        jetflavor = 1;
      }
      if(jettemp.DeltaR(truth_Z)<0.8){
        jetflavor = 2;
      }

      //running basic tjets
      tjetvar_1axis = tjet(1,inclusive_jets[ijet]);
      tjetvar_2axis = tjet(2,inclusive_jets[ijet]);

      Truth_flavor = jetflavor;
      Truth_pt     = jettemp.Pt();
      Truth_eta    = jettemp.Eta();
      Truth_phi    = jettemp.Phi();
      Truth_m      = jettemp.M();
      Truth_tjet1  = tjetvar_1axis;
      Truth_tjet2  = tjetvar_2axis;
      Truth_tau1   = nsub1(inclusive_jets[ijet]);
      Truth_tau2   = nsub2(inclusive_jets[ijet]);
      Truth_c2     = ecfC2(inclusive_jets[ijet]);
      Truth_d2     = ecfD2(inclusive_jets[ijet]);

      cout<<"FillingJet: flav="<<Truth_flavor<<"  pt="<<Truth_pt<<"  m="<<Truth_m<<"  tjet1="<<Truth_tjet1<<" tjet2= "<<Truth_tjet2<<endl;

      treeout->Fill();
    }

  }


  TString outfilename="ntuple_WprimeWZ.root";
  fileout = new TFile( outfilename ,"RECREATE");
  treeout->Write();
  fileout->Close();







  return 0;

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

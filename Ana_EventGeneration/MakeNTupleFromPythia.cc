/*----------------------------------------------------------------------

TITLE : MakeNTupleFromPythia.cc

DESCRIPTION : This runs the Pythia8 event generator, parsing the particle list
to determine which particles are to be used as final state particles [status = 81-99]
and writes these out as a vector of doubles for the components of the 4-vectors.  It further
writes out 4-vectors for each of the truth particles coming from the process itself (partons, W/Z, top)
depending on the process being simulated that can be used for truth identification of jets after clustering

COMPILE :
$ source compile.sh

RUN :
$ ./MakeNTupleFromPythia <type> <output>

type   : 0 = dijet , 1 = ttbar , 2 = W'->WZ (this determines which config file is used)
output : Anything you want - but being logical.  This will be the input for the NTupler code

//----------------------------------------------------------------------*/

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <vector>
#include <string>

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main(int argc, char* argv[]){

  //debug flag
  bool debug=false;

  //exit if you dont pass a run card
  if(argc<3){
    cout<<"You need to specify more arguments"<<endl;
    return 1;
  }

  //get processtype
  int ProcessType = atoi(argv[1]);

  //get config file
  string ConfigFile   = "config_DEFAULT.cmnd";
  if(ProcessType==0)
    ConfigFile="config_dijet.cmnd";
  else if(ProcessType==1)
    ConfigFile="config_WprimeWZ.cmnd";
  else if(ProcessType==2)
    ConfigFile="config_ttbar.cmnd";
  else{
    cout<<"Bad process type!"<<endl;
    return 0;
  }

  //get output file name
  string OutfileName = argv[2];

  cout<<"InputArguments:  ProcessType="<<ProcessType<<"   ConfigFile="<<ConfigFile<<"  OutfileName="<<OutfileName<<endl;


  ///////////////////////////////////
  //initialize the output ttree and branches
  ///////////////////////////////////
  double truth_quark1_pt;
  double truth_quark1_eta;
  double truth_quark1_phi;
  double truth_quark1_m;

  double truth_quark2_pt;
  double truth_quark2_eta;
  double truth_quark2_phi;
  double truth_quark2_m;

  double truth_top1_pt;
  double truth_top1_eta;
  double truth_top1_phi;
  double truth_top1_m;

  double truth_top2_pt;
  double truth_top2_eta;
  double truth_top2_phi;
  double truth_top2_m;

  double truth_W_pt;
  double truth_W_eta;
  double truth_W_phi;
  double truth_W_m;

  double truth_Z_pt;
  double truth_Z_eta;
  double truth_Z_phi;
  double truth_Z_m;

  vector<int>    fspart_id;
  vector<double> fspart_pt;
  vector<double> fspart_eta;
  vector<double> fspart_phi;
  vector<double> fspart_m;

  TTree *tree = new TTree("tree","tree");
  tree->Branch("fspart_id", &fspart_id);
  tree->Branch("fspart_pt", &fspart_pt);
  tree->Branch("fspart_eta",&fspart_eta);
  tree->Branch("fspart_phi",&fspart_phi);
  tree->Branch("fspart_m",  &fspart_m);

  tree->Branch("truth_quark1_pt",  &truth_quark1_pt);
  tree->Branch("truth_quark1_eta", &truth_quark1_eta);
  tree->Branch("truth_quark1_phi", &truth_quark1_phi);
  tree->Branch("truth_quark1_m",   &truth_quark1_m);

  tree->Branch("truth_quark2_pt",  &truth_quark2_pt);
  tree->Branch("truth_quark2_eta", &truth_quark2_eta);
  tree->Branch("truth_quark2_phi", &truth_quark2_phi);
  tree->Branch("truth_quark2_m",   &truth_quark2_m);

  tree->Branch("truth_top1_pt",  &truth_top1_pt);
  tree->Branch("truth_top1_eta", &truth_top1_eta);
  tree->Branch("truth_top1_phi", &truth_top1_phi);
  tree->Branch("truth_top1_m",   &truth_top1_m);

  tree->Branch("truth_top2_pt",  &truth_top2_pt);
  tree->Branch("truth_top2_eta", &truth_top2_eta);
  tree->Branch("truth_top2_phi", &truth_top2_phi);
  tree->Branch("truth_top2_m",   &truth_top2_m);

  tree->Branch("truth_W_pt",  &truth_W_pt);
  tree->Branch("truth_W_eta", &truth_W_eta);
  tree->Branch("truth_W_phi", &truth_W_phi);
  tree->Branch("truth_W_m",   &truth_W_m);

  tree->Branch("truth_Z_pt",  &truth_Z_pt);
  tree->Branch("truth_Z_eta", &truth_Z_eta);
  tree->Branch("truth_Z_phi", &truth_Z_phi);
  tree->Branch("truth_Z_m",   &truth_Z_m);

  //initialize pythia
  Pythia pythia;

  //read in config card
  pythia.readFile(ConfigFile);

  //get the number of events from the config file
  int nEvents = pythia.mode("Main:numberOfEvents");

  //collide protons with setup defined above
  pythia.init();

  int nAcceptedEvents=0;

  ///////////////////////////////////
  //main event loop
  ///////////////////////////////////
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    //print out to show progress
    if(iEvent%100==0)
      cout<<"Event: "<<iEvent<<endl;

    //generates a new event
    pythia.next();

    //accept event flag
    bool acceptevent=false;

    //clear all branches in preparation for next event
    fspart_id.clear();
    fspart_pt.clear();
    fspart_eta.clear();
    fspart_phi.clear();
    fspart_m.clear();

    truth_quark1_pt  = 100.0;
    truth_quark1_eta = 10.0;
    truth_quark1_phi = 0.0;
    truth_quark1_m   = 1.0;

    truth_quark2_pt  = 100.0;
    truth_quark2_eta = 10.0;
    truth_quark2_phi = 0.0;
    truth_quark2_m   = 1.0;

    truth_W_pt  = 100.0;
    truth_W_eta = 10.0;
    truth_W_phi = 0.0;
    truth_W_m   = 1.0;

    truth_Z_pt  = 100.0;
    truth_Z_eta = 10.0;
    truth_Z_phi = 0.0;
    truth_Z_m   = 1.0;

    truth_top1_pt  = 100.0;
    truth_top1_eta = 10.0;
    truth_top1_phi = 0.0;
    truth_top1_m   = 1.0;

    truth_top2_pt  = 100.0;
    truth_top2_eta = 10.0;
    truth_top2_phi = 0.0;
    truth_top2_m   = 1.0;

    //loops through the particles in the event just generated
    for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
      //event filter for
      if(debug) cout<<"Process specific filters: "<<ProcessType<<endl;
      if(ProcessType==0){
        if(debug) cout<<"Dijet Filter"<<endl;
        acceptevent=true;
        //fill truth quark branches
        if(pythia.event[iPart].id()<0 && pythia.event[iPart].status()==23){
          truth_quark1_pt=pythia.event[iPart].pT();
          truth_quark1_eta=pythia.event[iPart].eta();
          truth_quark1_phi=pythia.event[iPart].phi();
          truth_quark1_m=pythia.event[iPart].m();
        }
        if(pythia.event[iPart].id()>0 && pythia.event[iPart].status()==23){
          truth_quark2_pt=pythia.event[iPart].pT();
          truth_quark2_eta=pythia.event[iPart].eta();
          truth_quark2_phi=pythia.event[iPart].phi();
          truth_quark2_m=pythia.event[iPart].m();
        }
      }
      else if(ProcessType==1){
        if(debug) cout<<"WZ Filter"<<endl;
        //only accept if Wprime decay is via W'->WZ - identify by asking for Z0 in the intermediate state
        if(pythia.event[iPart].id()==23 && pythia.event[iPart].status()==-22){
          acceptevent=true;
        }
        //fill truth boson branches
        if(pythia.event[iPart].id()==23 && pythia.event[iPart].status()==-22){
          truth_Z_pt  = pythia.event[iPart].pT();
          truth_Z_eta = pythia.event[iPart].eta();
          truth_Z_phi = pythia.event[iPart].phi();
          truth_Z_m   = pythia.event[iPart].m();
        }
        if( (pythia.event[iPart].id()==24 || pythia.event[iPart].id()==-24) && pythia.event[iPart].status()==-22){
          truth_W_pt  = pythia.event[iPart].pT();
          truth_W_eta = pythia.event[iPart].eta();
          truth_W_phi = pythia.event[iPart].phi();
          truth_W_m   = pythia.event[iPart].m();
        }
      }
      else if(ProcessType=2){
        if(debug) cout<<"Top quark Filter"<<endl;
        acceptevent=true;
        //fill truth topquark branches
        if(pythia.event[iPart].id()==-6 && pythia.event[iPart].status()==-62){
          truth_top1_pt=pythia.event[iPart].pT();
          truth_top1_eta=pythia.event[iPart].eta();
          truth_top1_phi=pythia.event[iPart].phi();
          truth_top1_m=pythia.event[iPart].m();
        }
        if(pythia.event[iPart].id()==6 && pythia.event[iPart].status()==-62){
          truth_top2_pt=pythia.event[iPart].pT();
          truth_top2_eta=pythia.event[iPart].eta();
          truth_top2_phi=pythia.event[iPart].phi();
          truth_top2_m=pythia.event[iPart].m();
        }
      }
      else{
        acceptevent=true;
      }
      if(debug) cout<<acceptevent<<endl;

      //only save particles to output ttree if they are final state particles
      if(pythia.event[iPart].status()>=81 && pythia.event[iPart].status()<=99){

        if(debug) {
          if(iEvent%10==0){
            cout<<pythia.event[iPart].status()<<"  "
                <<pythia.event[iPart].id()<<"  "
                <<pythia.event[iPart].pT()<<"  "
                <<pythia.event[iPart].eta()<<"  "
                <<pythia.event[iPart].phi()<<"  "
                <<pythia.event[iPart].m()<<endl;
          }
        }

        fspart_id .push_back(pythia.event[iPart].id());
        fspart_pt .push_back(pythia.event[iPart].pT());
        fspart_eta.push_back(pythia.event[iPart].eta());
        fspart_phi.push_back(pythia.event[iPart].phi());
        fspart_m  .push_back(pythia.event[iPart].m());
      }
    }

    //fill event into tree if it has passed the event filter
    if(acceptevent){
      nAcceptedEvents+=1;
      tree->Fill();
    }


  }

  ///////////////////////////////////
  //Summary of truth filters
  ///////////////////////////////////
  cout<<"ProcessedEvents:"<<endl
      <<"total      = "<<nEvents<<endl
      <<"accepted   = "<<nAcceptedEvents<<endl
      <<"efficiency = "<<float(nAcceptedEvents)/float(nEvents)<<endl;

  TFile *fout = new TFile(OutfileName.c_str(),"RECREATE");
  tree->Write("tree");
  fout->Close();

  return 0;
}

#ifndef NTUPLER_H
#define NTUPLER_H

#define CERRD cout<<"Problem on "<<__FILE__<<"  "<<__LINE__<<endl;

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>

#include "TelescopingJets.hh"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

using namespace fastjet;
using namespace std;


///////////////////////////
//input file and tree
///////////////////////////
TTree *treein;
TFile *filein;

int nEvents;

///////////////////////////
//input tree branches
///////////////////////////
double truth_q1_pt;
double truth_q1_eta;
double truth_q1_phi;
double truth_q1_m;

double truth_q2_pt;
double truth_q2_eta;
double truth_q2_phi;
double truth_q2_m;

double truth_t1_pt;
double truth_t1_eta;
double truth_t1_phi;
double truth_t1_m;

double truth_t2_pt;
double truth_t2_eta;
double truth_t2_phi;
double truth_t2_m;

double truth_W_pt;
double truth_W_eta;
double truth_W_phi;
double truth_W_m;

double truth_Z_pt;
double truth_Z_eta;
double truth_Z_phi;
double truth_Z_m;

vector<int>*    fspart_id;
vector<double>* fspart_pt;
vector<double>* fspart_eta;
vector<double>* fspart_phi;
vector<double>* fspart_m;



///////////////////////////
//output file and tree
///////////////////////////
TTree *treeout;
TFile *fileout;

///////////////////////////
//output tree branches
///////////////////////////
int    NumberOfVertices;

int    TruthRaw_flavor;
double TruthRaw_pt;
double TruthRaw_eta;
double TruthRaw_phi;
double TruthRaw_m;
double TruthRaw_Tau1;
double TruthRaw_Tau2;
double TruthRaw_Tau21;
double TruthRaw_C2;
double TruthRaw_D2;
double TruthRaw_C3;
double TruthRaw_TJet_m1;
double TruthRaw_TJet_m2;
double TruthRaw_TJet_Tau1;
double TruthRaw_TJet_Tau2;
double TruthRaw_TJet_Tau21;
double TruthRaw_TJet_C2;
double TruthRaw_TJet_D2;
double TruthRaw_TJet_C3;

int    TruthPileup_flavor;
double TruthPileup_pt;
double TruthPileup_eta;
double TruthPileup_phi;
double TruthPileup_m;
double TruthPileup_Tau1;
double TruthPileup_Tau2;
double TruthPileup_Tau21;
double TruthPileup_C2;
double TruthPileup_D2;
double TruthPileup_C3;
double TruthPileup_TJet_m1;
double TruthPileup_TJet_m2;
double TruthPileup_TJet_Tau1;
double TruthPileup_TJet_Tau2;
double TruthPileup_TJet_Tau21;
double TruthPileup_TJet_C2;
double TruthPileup_TJet_D2;
double TruthPileup_TJet_C3;

int    RecoRaw_flavor;
double RecoRaw_pt;
double RecoRaw_eta;
double RecoRaw_phi;
double RecoRaw_m;
double RecoRaw_Tau1;
double RecoRaw_Tau2;
double RecoRaw_Tau21;
double RecoRaw_C2;
double RecoRaw_D2;
double RecoRaw_C3;
double RecoRaw_TJet_m1;
double RecoRaw_TJet_m2;
double RecoRaw_TJet_Tau1;
double RecoRaw_TJet_Tau2;
double RecoRaw_TJet_Tau21;
double RecoRaw_TJet_C2;
double RecoRaw_TJet_D2;
double RecoRaw_TJet_C3;

int    RecoPileup_flavor;
double RecoPileup_pt;
double RecoPileup_eta;
double RecoPileup_phi;
double RecoPileup_m;
double RecoPileup_Tau1;
double RecoPileup_Tau2;
double RecoPileup_Tau21;
double RecoPileup_C2;
double RecoPileup_D2;
double RecoPileup_C3;
double RecoPileup_TJet_m1;
double RecoPileup_TJet_m2;
double RecoPileup_TJet_Tau1;
double RecoPileup_TJet_Tau2;
double RecoPileup_TJet_Tau21;
double RecoPileup_TJet_C2;
double RecoPileup_TJet_D2;
double RecoPileup_TJet_C3;

///////////////////////////
//extra functions
///////////////////////////
vector<PseudoJet> ToyCalorimeter(vector<PseudoJet> truth_particles);
double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max);
double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max);

#endif

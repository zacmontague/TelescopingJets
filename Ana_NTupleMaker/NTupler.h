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
//
///////////////////////////
int    tempJet_flavor;
double tempJet_pt;
double tempJet_eta;
double tempJet_phi;
double tempJet_m;
double tempJet_Tau1;
double tempJet_Tau2;
double tempJet_Tau21;
double tempJet_C2;
double tempJet_D2;
double tempJet_C3;
double tempJet_TJet_m1;
double tempJet_TJet_m2;
double tempJet_TJet_Tau1;
double tempJet_TJet_Tau2;
double tempJet_TJet_Tau21;
double tempJet_TJet_C2;
double tempJet_TJet_D2;
double tempJet_TJet_C3;

///////////////////////////
//output tree branches
///////////////////////////
int    NumberOfVertices;

vector<int>    TruthRaw_flavor;
vector<double> TruthRaw_pt;
vector<double> TruthRaw_eta;
vector<double> TruthRaw_phi;
vector<double> TruthRaw_m;
vector<double> TruthRaw_Tau1;
vector<double> TruthRaw_Tau2;
vector<double> TruthRaw_Tau21;
vector<double> TruthRaw_C2;
vector<double> TruthRaw_D2;
vector<double> TruthRaw_C3;
vector<double> TruthRaw_TJet_m1;
vector<double> TruthRaw_TJet_m2;
vector<double> TruthRaw_TJet_Tau1;
vector<double> TruthRaw_TJet_Tau2;
vector<double> TruthRaw_TJet_Tau21;
vector<double> TruthRaw_TJet_C2;
vector<double> TruthRaw_TJet_D2;
vector<double> TruthRaw_TJet_C3;

vector<int>    TruthPileup_flavor;
vector<double> TruthPileup_pt;
vector<double> TruthPileup_eta;
vector<double> TruthPileup_phi;
vector<double> TruthPileup_m;
vector<double> TruthPileup_Tau1;
vector<double> TruthPileup_Tau2;
vector<double> TruthPileup_Tau21;
vector<double> TruthPileup_C2;
vector<double> TruthPileup_D2;
vector<double> TruthPileup_C3;
vector<double> TruthPileup_TJet_m1;
vector<double> TruthPileup_TJet_m2;
vector<double> TruthPileup_TJet_Tau1;
vector<double> TruthPileup_TJet_Tau2;
vector<double> TruthPileup_TJet_Tau21;
vector<double> TruthPileup_TJet_C2;
vector<double> TruthPileup_TJet_D2;
vector<double> TruthPileup_TJet_C3;

vector<int>    RecoRaw_flavor;
vector<double> RecoRaw_pt;
vector<double> RecoRaw_eta;
vector<double> RecoRaw_phi;
vector<double> RecoRaw_m;
vector<double> RecoRaw_Tau1;
vector<double> RecoRaw_Tau2;
vector<double> RecoRaw_Tau21;
vector<double> RecoRaw_C2;
vector<double> RecoRaw_D2;
vector<double> RecoRaw_C3;
vector<double> RecoRaw_TJet_m1;
vector<double> RecoRaw_TJet_m2;
vector<double> RecoRaw_TJet_Tau1;
vector<double> RecoRaw_TJet_Tau2;
vector<double> RecoRaw_TJet_Tau21;
vector<double> RecoRaw_TJet_C2;
vector<double> RecoRaw_TJet_D2;
vector<double> RecoRaw_TJet_C3;

vector<int>    RecoPileup_flavor;
vector<double> RecoPileup_pt;
vector<double> RecoPileup_eta;
vector<double> RecoPileup_phi;
vector<double> RecoPileup_m;
vector<double> RecoPileup_Tau1;
vector<double> RecoPileup_Tau2;
vector<double> RecoPileup_Tau21;
vector<double> RecoPileup_C2;
vector<double> RecoPileup_D2;
vector<double> RecoPileup_C3;
vector<double> RecoPileup_TJet_m1;
vector<double> RecoPileup_TJet_m2;
vector<double> RecoPileup_TJet_Tau1;
vector<double> RecoPileup_TJet_Tau2;
vector<double> RecoPileup_TJet_Tau21;
vector<double> RecoPileup_TJet_C2;
vector<double> RecoPileup_TJet_D2;
vector<double> RecoPileup_TJet_C3;

vector<int>    TruthRawTrim_flavor;
vector<double> TruthRawTrim_pt;
vector<double> TruthRawTrim_eta;
vector<double> TruthRawTrim_phi;
vector<double> TruthRawTrim_m;
vector<double> TruthRawTrim_Tau1;
vector<double> TruthRawTrim_Tau2;
vector<double> TruthRawTrim_Tau21;
vector<double> TruthRawTrim_C2;
vector<double> TruthRawTrim_D2;
vector<double> TruthRawTrim_C3;
vector<double> TruthRawTrim_TJet_m1;
vector<double> TruthRawTrim_TJet_m2;
vector<double> TruthRawTrim_TJet_Tau1;
vector<double> TruthRawTrim_TJet_Tau2;
vector<double> TruthRawTrim_TJet_Tau21;
vector<double> TruthRawTrim_TJet_C2;
vector<double> TruthRawTrim_TJet_D2;
vector<double> TruthRawTrim_TJet_C3;

vector<int>    TruthPileupTrim_flavor;
vector<double> TruthPileupTrim_pt;
vector<double> TruthPileupTrim_eta;
vector<double> TruthPileupTrim_phi;
vector<double> TruthPileupTrim_m;
vector<double> TruthPileupTrim_Tau1;
vector<double> TruthPileupTrim_Tau2;
vector<double> TruthPileupTrim_Tau21;
vector<double> TruthPileupTrim_C2;
vector<double> TruthPileupTrim_D2;
vector<double> TruthPileupTrim_C3;
vector<double> TruthPileupTrim_TJet_m1;
vector<double> TruthPileupTrim_TJet_m2;
vector<double> TruthPileupTrim_TJet_Tau1;
vector<double> TruthPileupTrim_TJet_Tau2;
vector<double> TruthPileupTrim_TJet_Tau21;
vector<double> TruthPileupTrim_TJet_C2;
vector<double> TruthPileupTrim_TJet_D2;
vector<double> TruthPileupTrim_TJet_C3;

vector<int>    RecoRawTrim_flavor;
vector<double> RecoRawTrim_pt;
vector<double> RecoRawTrim_eta;
vector<double> RecoRawTrim_phi;
vector<double> RecoRawTrim_m;
vector<double> RecoRawTrim_Tau1;
vector<double> RecoRawTrim_Tau2;
vector<double> RecoRawTrim_Tau21;
vector<double> RecoRawTrim_C2;
vector<double> RecoRawTrim_D2;
vector<double> RecoRawTrim_C3;
vector<double> RecoRawTrim_TJet_m1;
vector<double> RecoRawTrim_TJet_m2;
vector<double> RecoRawTrim_TJet_Tau1;
vector<double> RecoRawTrim_TJet_Tau2;
vector<double> RecoRawTrim_TJet_Tau21;
vector<double> RecoRawTrim_TJet_C2;
vector<double> RecoRawTrim_TJet_D2;
vector<double> RecoRawTrim_TJet_C3;

vector<int>    RecoPileupTrim_flavor;
vector<double> RecoPileupTrim_pt;
vector<double> RecoPileupTrim_eta;
vector<double> RecoPileupTrim_phi;
vector<double> RecoPileupTrim_m;
vector<double> RecoPileupTrim_Tau1;
vector<double> RecoPileupTrim_Tau2;
vector<double> RecoPileupTrim_Tau21;
vector<double> RecoPileupTrim_C2;
vector<double> RecoPileupTrim_D2;
vector<double> RecoPileupTrim_C3;
vector<double> RecoPileupTrim_TJet_m1;
vector<double> RecoPileupTrim_TJet_m2;
vector<double> RecoPileupTrim_TJet_Tau1;
vector<double> RecoPileupTrim_TJet_Tau2;
vector<double> RecoPileupTrim_TJet_Tau21;
vector<double> RecoPileupTrim_TJet_C2;
vector<double> RecoPileupTrim_TJet_D2;
vector<double> RecoPileupTrim_TJet_C3;

///////////////////////////
//extra functions
///////////////////////////
void ResetBranches();

vector<PseudoJet> ToyCalorimeter(vector<PseudoJet> truth_particles);
double T_Nsubjettiness(int N, PseudoJet& input, double beta_min, double beta_max);
double T_NsubjettinessRatio(int N_num, int N_den, PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_C2(PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_D2(PseudoJet& input, double beta_min, double beta_max);
double T_EnergyCorrelator_C3(PseudoJet& input, double beta_min, double beta_max);

#endif

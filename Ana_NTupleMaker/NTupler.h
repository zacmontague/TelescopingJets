#ifndef NTUPLER_H
#define NTUPLER_H

#define CERRD cout<<"Problem on "<<__FILE__<<"  "<<__LINE__<<endl;

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TLorentzVector.h"

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





using namespace fastjet;
using namespace std;


//input file and tree
TTree *treein;
TFile *filein;

int nEvents;

///////////////////////////
//input tree branches
///////////////////////////
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

vector<int>*    fspart_id;
vector<double>* fspart_pt;
vector<double>* fspart_eta;
vector<double>* fspart_phi;
vector<double>* fspart_m;



//temp variables
double tjetvar_1axis;
double tjetvar_2axis;
int jetvariable;

//output file and tree
TTree *treeout;
TFile *fileout;

///////////////////////////
//output tree branches
///////////////////////////
int    Truth_flavor;
double Truth_pt;
double Truth_eta;
double Truth_phi;
double Truth_m;
double Truth_tjet1;
double Truth_tjet2;



#endif
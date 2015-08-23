#include <TCanvas.h>

void RunAllVerify(){

  TFile *f0 = new TFile("Data/ntuple_dijet_single.root");
  TFile *f1 = new TFile("Data/ntuple_tt_single.root");
  TFile *f2 = new TFile("Data/ntuple_ww_single.root");
  TFile *f3 = new TFile("Data/ntuple_zz_single.root");
  TFile *f4 = new TFile("Data/ntuple_hh_single.root");

  TTree *t0 = f0->Get("JetTree");
  TTree *t1 = f1->Get("JetTree");
  TTree *t2 = f2->Get("JetTree");
  TTree *t3 = f3->Get("JetTree");
  TTree *t4 = f4->Get("JetTree");


  ///////////////////////////
  //SANITY CHECK
  ///////////////////////////

  t0->Draw("TruthRaw_m>>h0(60,0,300)");
  t1->Draw("TruthRaw_m>>h1(60,0,300)");
  t2->Draw("TruthRaw_m>>h2(60,0,300)");
  t3->Draw("TruthRaw_m>>h3(60,0,300)");
  t4->Draw("TruthRaw_m>>h4(60,0,300)");

  TH1D *h0 = h0->Clone("h0");
  TH1D *h1 = h1->Clone("h1");
  TH1D *h2 = h2->Clone("h2");
  TH1D *h3 = h3->Clone("h3");
  TH1D *h4 = h4->Clone("h4");

  h0->Scale(1.0/h0->Integral());
  h1->Scale(1.0/h1->Integral());
  h2->Scale(1.0/h2->Integral());
  h3->Scale(1.0/h3->Integral());
  h4->Scale(1.0/h4->Integral());

  h0->SetMaximum(0.3);
  h1->SetMaximum(0.3);
  h2->SetMaximum(0.3);
  h3->SetMaximum(0.3);
  h4->SetMaximum(0.3);

  h0->GetXaxis()->SetTitle("M(jet) [GeV]");
  h1->GetXaxis()->SetTitle("M(jet) [GeV]");
  h2->GetXaxis()->SetTitle("M(jet) [GeV]");
  h3->GetXaxis()->SetTitle("M(jet) [GeV]");
  h4->GetXaxis()->SetTitle("M(jet) [GeV]");

  h0->GetYaxis()->SetTitle("Normalized Units");
  h1->GetYaxis()->SetTitle("Normalized Units");
  h2->GetYaxis()->SetTitle("Normalized Units");
  h3->GetYaxis()->SetTitle("Normalized Units");
  h4->GetYaxis()->SetTitle("Normalized Units");

  h0->SetLineColor(1);
  h1->SetLineColor(2);
  h2->SetLineColor(3);
  h3->SetLineColor(4);
  h4->SetLineColor(9);

  h0->SetTitle("Dijet(qq,qg,gg)");
  h1->SetTitle("G#rightarrow tt#rightarrow WbWb#rightarrow qqbqqb");
  h2->SetTitle("G#rightarrow WW#rightarrow qqqq");
  h3->SetTitle("G#rightarrow ZZ#rightarrow qqqq");
  h4->SetTitle("G#rightarrow HH#rightarrow bbbb");

  TCanvas *can = new TCanvas("can","can",400,800);
  can->SetTitle("Sanity Check");
  can->Divide(1,5);
  can->cd(1);
  h0->Draw();
  can->cd(2);
  h1->Draw();
  can->cd(3);
  h2->Draw();
  can->cd(4);
  h3->Draw();
  can->cd(5);
  h4->Draw();


  ///////////////////////////
  //TRUTH MATCHING CHECK
  ///////////////////////////
  t1->Draw("TruthRaw_m>>ht_notmatched(60,0,300)","TruthRaw_flavor!=3");
  t1->Draw("TruthRaw_m>>ht_matched(60,0,300)","TruthRaw_flavor==3");
  t2->Draw("TruthRaw_m>>hw_notmatched(60,0,300)","TruthRaw_flavor!=1");
  t2->Draw("TruthRaw_m>>hw_matched(60,0,300)","TruthRaw_flavor==1");
  t3->Draw("TruthRaw_m>>hz_notmatched(60,0,300)","TruthRaw_flavor!=2");
  t3->Draw("TruthRaw_m>>hz_matched(60,0,300)","TruthRaw_flavor==2");
  t4->Draw("TruthRaw_m>>hh_notmatched(60,0,300)","TruthRaw_flavor!=4");
  t4->Draw("TruthRaw_m>>hh_matched(60,0,300)","TruthRaw_flavor==4");

  TH1D* ht_notmatched = ht_notmatched->Clone("ht_notmatched");
  TH1D* ht_matched = ht_matched->Clone("ht_matched");
  TH1D* hw_notmatched = hw_notmatched->Clone("hw_notmatched");
  TH1D* hw_matched = hw_matched->Clone("hw_matched");
  TH1D* hz_notmatched = hz_notmatched->Clone("hz_notmatched");
  TH1D* hz_matched = hz_matched->Clone("hz_matched");
  TH1D* hh_notmatched = hh_notmatched->Clone("hh_notmatched");
  TH1D* hh_matched = hh_matched->Clone("hh_matched");

  ht_notmatched->SetLineColor(2);
  ht_matched   ->SetLineColor(2);
  hw_notmatched->SetLineColor(3);
  hw_matched   ->SetLineColor(3);
  hz_notmatched->SetLineColor(4);
  hz_matched   ->SetLineColor(4);
  hh_notmatched->SetLineColor(9);
  hh_matched   ->SetLineColor(9);

  ht_notmatched->SetTitle("ttbar NOT Matched");
  ht_matched   ->SetTitle("ttbar Matched");
  hw_notmatched->SetTitle("W NOT Matched");
  hw_matched   ->SetTitle("W Matched");
  hz_notmatched->SetTitle("Z NOT Matched");
  hz_matched   ->SetTitle("Z Matched");
  hh_notmatched->SetTitle("h NOT Matched");
  hh_matched   ->SetTitle("h Matched");

  TCanvas *can1 = new TCanvas("can1","can1",1000,400);
  can1->SetTitle("Truth Matching Check");
  can1->Divide(4,2);
  can1->cd(1);
  ht_notmatched->Draw();
  can1->cd(5);
  ht_matched->Draw();

  can1->cd(2);
  hw_notmatched->Draw();
  can1->cd(6);
  hw_matched->Draw();

  can1->cd(3);
  hz_notmatched->Draw();
  can1->cd(7);
  hz_matched->Draw();

  can1->cd(4);
  hh_notmatched->Draw();
  can1->cd(8);
  hh_matched->Draw();



}

#include <TROOT.h>
#include <TPad.h>
#include <TFrame.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <cmath>

void makeplot (string isoregion, string plotname) {

  //definition of stacking plotname
  THStack *hs = new THStack("hs","Stacked histograms");

  //opening background files
  TFile *allDY_ifile = TFile::Open("crvalidation/allDY.root");
  TFile *allQCD_ifile = TFile::Open("crvalidation/allQCD.root");
  TFile *allTT_ifile = TFile::Open("crvalidation/allTT.root");
  TFile *allVV_ifile = TFile::Open("crvalidation/allVV.root");
  TFile *allWJets_ifile = TFile::Open("crvalidation/allWJets.root");

  //defining the different background plots
  TH1F* crplot_allDY;
  TH1F* crplot_allQCD;
  TH1F* crplot_allTT;
  TH1F* crplot_allVV;
  TH1F* crplot_allWJets;
  TH1F* crplot_allSum;
  TH1F* crplot_allSum_err;

  crplot_allDY = ((TH1F*)(allDY_ifile->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_allDY->SetFillColor(594);
  crplot_allDY->SetLineColor(594);
  crplot_allDY->SetStats(0);

  crplot_allQCD = ((TH1F*)(allQCD_ifile->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_allQCD->SetFillColor(426);
  crplot_allQCD->SetLineColor(426);
  crplot_allQCD->SetStats(0);
  crplot_allQCD->GetYaxis()->SetRangeUser(0.,10000000.);

  crplot_allTT = ((TH1F*)(allTT_ifile->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_allTT->SetFillColor(401);
  crplot_allTT->SetLineColor(401);
  crplot_allTT->SetStats(0);

  crplot_allVV = ((TH1F*)(allVV_ifile->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_allVV->SetFillColor(626);
  crplot_allVV->SetLineColor(626);
  crplot_allVV->SetStats(0);

  crplot_allWJets = ((TH1F*)(allWJets_ifile->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_allWJets->SetFillColor(610);
  crplot_allWJets->SetLineColor(610);
  crplot_allWJets->SetStats(0);

  //creating and rescaling the bkg sum plot
  crplot_allSum = (TH1F*) crplot_allDY->Clone();
  crplot_allSum->SetTitle("CMS Work in Progress");
  crplot_allSum->Add(crplot_allQCD);
  crplot_allSum->Add(crplot_allTT);
  crplot_allSum->Add(crplot_allVV);
  crplot_allSum->Add(crplot_allWJets);

  double maximum = 0.;
  double minimum = FLT_MAX;
  int nbins = crplot_allSum->GetNbinsX();

	for (int i = 0; i < nbins; i++) {
    if ( (crplot_allSum->GetBinContent(i) < minimum) && (crplot_allSum->GetBinContent(i) != 0.) ) minimum = crplot_allSum->GetBinContent(i);
    if ( crplot_allSum->GetBinContent(i) > maximum) maximum = crplot_allSum->GetBinContent(i);
  }

  crplot_allSum->GetYaxis()->SetRangeUser( (minimum * 0.1) , (maximum * 10.) );

  //extra plot cosmetics
  if ( plotname == "h_dijetinvariantmass" ) {
    crplot_allSum->GetXaxis()->SetTitle("M_{(jet,jet)} [GeV]");
  }

  //defining errorbar plot
  crplot_allSum_err = (TH1F*) crplot_allSum->Clone();
  crplot_allSum_err->SetFillColor(kRed);
  crplot_allSum_err->SetFillStyle(3002);
  crplot_allSum_err->SetLineColor(kRed);


  //defining legend
	TLegend* leg = new TLegend(0.70,0.65,0.9,0.9);
  leg->AddEntry(crplot_allDY,"Drell-Yan","f");
  leg->AddEntry(crplot_allQCD,"QCD","f");
  leg->AddEntry(crplot_allTT,"t#bar{t}","f");
  leg->AddEntry(crplot_allVV,"VV","f");
  leg->AddEntry(crplot_allWJets,"W + jets","f");
  leg->AddEntry(crplot_allSum_err, "#sigma^{total}_{stat}");

  //creating canvas
  TCanvas *my_canvas = new TCanvas("c1","c1",10,10,600,600);
  //TCanvas *my_canvas = new TCanvas;
	my_canvas->cd();

  //stacking plots
  hs->Add(crplot_allDY);
  hs->Add(crplot_allVV);
  hs->Add(crplot_allTT);
  hs->Add(crplot_allWJets);
  hs->Add(crplot_allQCD);

  //painting and saving plots into canvas
  crplot_allSum->Draw("HIST");
  hs->Draw("HIST+same");
  crplot_allSum_err->Draw("e2+same");
  leg->Draw();
  gPad->SetLogy();
  my_canvas->Print(("results/" + plotname + "_" + isoregion +".pdf").c_str());
  my_canvas->Close();
}

void printplots() {
  makeplot("Tau2TightIsoVBFInverted","h_dijetinvariantmass");
}

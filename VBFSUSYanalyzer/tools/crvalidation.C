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
  TFile *signal_lsp000_chi100 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root");
  TFile *signal_lsp000_chi200 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root");
  TFile *signal_lsp000_chi300 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
  TFile *signal_lsp050_chi100 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau095_Chargino100_1M.root");
  TFile *signal_lsp050_chi200 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau195_Chargino200_1M.root");
  TFile *signal_lsp050_chi300 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau295_Chargino300_1M.root");

  //defining the different background plots
  TH1F* crplot_allDY;
  TH1F* crplot_allQCD;
  TH1F* crplot_allTT;
  TH1F* crplot_allVV;
  TH1F* crplot_allWJets;
  TH1F* crplot_allSum;
  TH1F* crplot_allSum_err;

  TH1F* crplot_sig_lsp000_chi100 = ((TH1F*)(signal_lsp000_chi100->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp000_chi100->SetLineColor(kBlue);
  crplot_sig_lsp000_chi100->SetLineWidth(3);
  crplot_sig_lsp000_chi100->SetStats(0);

  TH1F* crplot_sig_lsp000_chi200 = ((TH1F*)(signal_lsp000_chi200->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp000_chi200->SetLineColor(kRed);
  crplot_sig_lsp000_chi200->SetLineWidth(3);
  crplot_sig_lsp000_chi200->SetStats(0);

  TH1F* crplot_sig_lsp000_chi300 = ((TH1F*)(signal_lsp000_chi300->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp000_chi300->SetLineColor(kBlack);
  crplot_sig_lsp000_chi300->SetLineWidth(3);
  crplot_sig_lsp000_chi300->SetStats(0);

  TH1F* crplot_sig_lsp050_chi100 = ((TH1F*)(signal_lsp050_chi100->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp050_chi100->SetLineColor(kBlue);
  crplot_sig_lsp050_chi100->SetLineWidth(3);
  crplot_sig_lsp050_chi100->SetStats(0);

  TH1F* crplot_sig_lsp050_chi200 = ((TH1F*)(signal_lsp050_chi200->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp050_chi200->SetLineColor(kRed);
  crplot_sig_lsp050_chi200->SetLineWidth(3);
  crplot_sig_lsp050_chi200->SetStats(0);

  TH1F* crplot_sig_lsp050_chi300 = ((TH1F*)(signal_lsp050_chi300->Get(("demo/" + isoregion +"ObjectSelection/" + plotname).c_str())));
  crplot_sig_lsp050_chi300->SetLineColor(kBlack);
  crplot_sig_lsp050_chi300->SetLineWidth(3);
  crplot_sig_lsp050_chi300->SetStats(0);

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

  //creating and Yaxis rescaling the bkg sum and sig plots
  crplot_allSum = (TH1F*) crplot_allDY->Clone();
  crplot_allSum->SetTitle("Simulation 13 TeV");
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

    if ( (crplot_sig_lsp000_chi100->GetBinContent(i) < minimum) && (crplot_sig_lsp000_chi100->GetBinContent(i) != 0.) ) minimum = crplot_sig_lsp000_chi100->GetBinContent(i);
    if ( crplot_sig_lsp000_chi100->GetBinContent(i) > maximum) maximum = crplot_sig_lsp000_chi100->GetBinContent(i);

    if ( (crplot_sig_lsp000_chi200->GetBinContent(i) < minimum) && (crplot_sig_lsp000_chi200->GetBinContent(i) != 0.) ) minimum = crplot_sig_lsp000_chi200->GetBinContent(i);
    if ( crplot_sig_lsp000_chi200->GetBinContent(i) > maximum) maximum = crplot_sig_lsp000_chi200->GetBinContent(i);

    if ( (crplot_sig_lsp000_chi300->GetBinContent(i) < minimum) && (crplot_sig_lsp000_chi300->GetBinContent(i) != 0.) ) minimum = crplot_sig_lsp000_chi300->GetBinContent(i);
    if ( crplot_sig_lsp000_chi300->GetBinContent(i) > maximum) maximum = crplot_sig_lsp000_chi300->GetBinContent(i);

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
	TLegend* leg = new TLegend(0.70,0.55,0.9,0.9);
  leg->SetTextSize(0.02);
  leg->AddEntry(crplot_allDY,"Drell-Yan","f");
  leg->AddEntry(crplot_allQCD,"QCD","f");
  leg->AddEntry(crplot_allTT,"t#bar{t}","f");
  leg->AddEntry(crplot_allVV,"VV","f");
  leg->AddEntry(crplot_allWJets,"W + jets","f");
  leg->AddEntry(crplot_allSum_err, "#sigma^{total}_{stat}");
  leg->AddEntry(crplot_sig_lsp000_chi100, "#tilde{#chi}^{#pm}_{1} = 100 GeV");
  leg->AddEntry(crplot_sig_lsp000_chi200, "#tilde{#chi}^{#pm}_{1} = 200 GeV");
  leg->AddEntry(crplot_sig_lsp000_chi300, "#tilde{#chi}^{#pm}_{1} = 300 GeV");
  leg->AddEntry((TObject*)0, "#tilde{#chi}^{0}_{1} = 0 GeV", "");

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
  crplot_sig_lsp000_chi100->Draw("HIST+same");
  crplot_sig_lsp000_chi200->Draw("HIST+same");
  crplot_sig_lsp000_chi300->Draw("HIST+same");
  leg->Draw();
  gPad->SetLogy();
  my_canvas->Print(("results/" + plotname + "_" + isoregion +".pdf").c_str());
  my_canvas->Close();
}

void printplots() {
  makeplot("Tau2TightIsoVBFInverted","h_dijetinvariantmass");
  makeplot("Tau2TightIsoVBFInverted","h_met");
  makeplot("Tau2TightIsoVBFInverted","h_tau2pt");
  makeplot("Tau2TightIsoVBFInverted","h_jet1pt");

  makeplot("Taui2TightIso","h_dijetinvariantmass");
  makeplot("Taui2TightIso","h_dijetdeltaeta");
  makeplot("Taui2TightIso","h_ditauinvariantmass");
  makeplot("Taui2TightIso","h_ditaudeltaeta");
  makeplot("Taui2TightIso","h_met");
  makeplot("Taui2TightIso","h_tau1pt");
  makeplot("Taui2TightIso","h_tau2pt");
  makeplot("Taui2TightIso","h_jet1pt");
  makeplot("Taui2TightIso","h_jet2pt");

  makeplot("Tau2LooseIsoInclusive","h_dijetinvariantmass");
  makeplot("Tau2LooseIsoInclusive","h_met");
  makeplot("Tau2LooseIsoInclusive","h_tau2pt");
  makeplot("Tau2LooseIsoInclusive","h_jet1pt");
  makeplot("Tau2LooseIsoInclusiveVBFInverted","h_dijetinvariantmass");
  makeplot("Tau2LooseIsoInclusiveVBFInverted","h_met");
  makeplot("Tau2LooseIsoInclusiveVBFInverted","h_tau2pt");
  makeplot("Tau2LooseIsoInclusiveVBFInverted","h_jet1pt");
}

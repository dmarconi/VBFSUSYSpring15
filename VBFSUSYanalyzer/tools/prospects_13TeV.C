#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <cmath>

void var_plot(string variable, string label, int rebin) {

	//loading input files
	TFile *inputfile_sig_1 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
	TFile *inputfile_sig_2 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root");
	TFile *inputfile_sig_3 = TFile::Open("crvalidation/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root");
	TFile *inputfile_qcd = TFile::Open("crvalidation/allQCD.root");
	TFile *inputfile_vv = TFile::Open("crvalidation/allVV.root");
	TFile *inputfile_wjets = TFile::Open("crvalidation/allWJets.root");

	//variables init
	Double_t norm = 1.;
	Double_t maximum = 0.;
	Double_t minimum = 10000.;
	Double_t temp_maximum = 0.;
	Double_t temp_minimum = 10000.;

	//definition and plot initialization
	TH1F* h_dijetinvariantmass_sig_1 = ((TH1F*)(inputfile_sig_1->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_sig_1->Scale(norm/h_dijetinvariantmass_sig_1->Integral());
	h_dijetinvariantmass_sig_1->SetTitle("Simulation, 13 TeV");
	h_dijetinvariantmass_sig_1->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_sig_1->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_sig_1->SetLineColor(kBlack);
	h_dijetinvariantmass_sig_1->SetLineWidth(3);
	h_dijetinvariantmass_sig_1->SetStats(kFALSE);
	h_dijetinvariantmass_sig_1->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_sig_1->GetMaximum();
	temp_minimum = h_dijetinvariantmass_sig_1->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_sig_2 = ((TH1F*)(inputfile_sig_2->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_sig_2->Scale(norm/h_dijetinvariantmass_sig_2->Integral());
	h_dijetinvariantmass_sig_2->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_sig_2->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_sig_2->SetLineColor(kRed);
	h_dijetinvariantmass_sig_2->SetLineWidth(3);
	h_dijetinvariantmass_sig_2->SetStats(kFALSE);
	h_dijetinvariantmass_sig_2->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_sig_2->GetMaximum();
	temp_minimum = h_dijetinvariantmass_sig_2->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_sig_3 = ((TH1F*)(inputfile_sig_3->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_sig_3->Scale(norm/h_dijetinvariantmass_sig_3->Integral());
	h_dijetinvariantmass_sig_3->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_sig_3->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_sig_3->SetLineColor(kBlue);
	h_dijetinvariantmass_sig_3->SetLineWidth(3);
	h_dijetinvariantmass_sig_3->SetStats(kFALSE);
	h_dijetinvariantmass_sig_3->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_sig_3->GetMaximum();
	temp_minimum = h_dijetinvariantmass_sig_3->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_qcd = ((TH1F*)(inputfile_qcd->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_qcd->SetLineColor(kCyan);
	h_dijetinvariantmass_qcd->SetFillColor(kCyan);
	h_dijetinvariantmass_qcd->SetFillStyle(3001);
	h_dijetinvariantmass_qcd->SetLineWidth(3);
	h_dijetinvariantmass_qcd->Scale(norm/h_dijetinvariantmass_qcd->Integral());
	h_dijetinvariantmass_qcd->SetTitle("");
	h_dijetinvariantmass_qcd->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_qcd->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_qcd->SetStats(kFALSE);
	h_dijetinvariantmass_qcd->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_qcd->GetMaximum();
	temp_minimum = h_dijetinvariantmass_qcd->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_vv = ((TH1F*)(inputfile_vv->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_vv->SetLineColor(kViolet+2);
	h_dijetinvariantmass_vv->SetFillColor(kViolet+2);
	h_dijetinvariantmass_vv->SetFillStyle(3001);
	h_dijetinvariantmass_vv->SetLineWidth(3);
	h_dijetinvariantmass_vv->Scale(norm/h_dijetinvariantmass_vv->Integral());
	h_dijetinvariantmass_vv->SetTitle("");
	h_dijetinvariantmass_vv->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_vv->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_vv->SetStats(kFALSE);
	h_dijetinvariantmass_vv->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_vv->GetMaximum();
	temp_minimum = h_dijetinvariantmass_vv->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;


	TH1F* h_dijetinvariantmass_wjets = ((TH1F*)(inputfile_wjets->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_wjets->SetLineColor(kGreen);
	h_dijetinvariantmass_wjets->SetFillColor(kGreen);
	h_dijetinvariantmass_wjets->SetFillStyle(3001);
	h_dijetinvariantmass_wjets->SetLineWidth(3);
	h_dijetinvariantmass_wjets->Scale(norm/h_dijetinvariantmass_wjets->Integral());
	h_dijetinvariantmass_wjets->SetTitle("");
	h_dijetinvariantmass_wjets->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_wjets->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_wjets->SetStats(kFALSE);
	h_dijetinvariantmass_wjets->Rebin(rebin);
	temp_maximum = h_dijetinvariantmass_wjets->GetMaximum();
	temp_minimum = h_dijetinvariantmass_wjets->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	//Setting proper Y axis range
	maximum = maximum + 0.3 * maximum;
	h_dijetinvariantmass_sig_1->GetYaxis()->SetRangeUser(minimum,maximum);
	h_dijetinvariantmass_qcd->GetYaxis()->SetRangeUser(minimum,maximum);
	h_dijetinvariantmass_vv->GetYaxis()->SetRangeUser(minimum,maximum);
	h_dijetinvariantmass_wjets->GetYaxis()->SetRangeUser(minimum,maximum);

	//defining legend
	TLegend* leg = new TLegend(0.58,0.7,0.9,0.9);
	leg->SetTextSize(0.02);
	leg->AddEntry(h_dijetinvariantmass_qcd,"QCD","f");
	leg->AddEntry(h_dijetinvariantmass_vv,"VVjj","f");
	leg->AddEntry(h_dijetinvariantmass_wjets,"W+jets","f");
	leg->AddEntry(h_dijetinvariantmass_sig_3,"#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}jj, #tilde{#chi}^{#pm}_{1} = 100 GeV","l");
	leg->AddEntry(h_dijetinvariantmass_sig_2,"#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}jj, #tilde{#chi}^{#pm}_{1} = 200 GeV","l");
	leg->AddEntry(h_dijetinvariantmass_sig_1,"#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}jj, #tilde{#chi}^{#pm}_{1} = 300 GeV","l");


	//plotting on canvas and writing to disk
	TCanvas *my_canvas = new TCanvas("my_canvas","my_canvas",500,500);
	my_canvas->cd();

	gPad->SetLogy();


	h_dijetinvariantmass_vv->Draw("HIST");
	h_dijetinvariantmass_wjets->Draw("HIST same");
	h_dijetinvariantmass_qcd->Draw("HIST same");
	h_dijetinvariantmass_sig_1->Draw("HIST same");
	h_dijetinvariantmass_sig_2->Draw("HIST same");
	h_dijetinvariantmass_sig_3->Draw("HIST same");
	leg->Draw();

	my_canvas->Print(("results/" + variable + "_prospects13tev.pdf").c_str());
	my_canvas->Close();

}

//main function
void prospects_13Tev (){
	var_plot("h_dijetinvariantmass", "m_{j_{1}j_{2}} [GeV]", 1);
	var_plot("h_tau1pt", "p_{t}(#tau_{1}) [GeV]", 5);
	var_plot("h_tau2pt", "p_{t}(#tau_{2}) [GeV]", 5);
	var_plot("h_jet1pt", "p_{t}(jet_{1}) [GeV]", 5);
	var_plot("h_met", "#slash{E}_{T} [GeV]", 1);
}

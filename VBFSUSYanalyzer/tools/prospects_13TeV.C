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

void var_plot(string variable, string label) {
	//loading input
	TFile *inputfile_sig = TFile::Open("prospects/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
	TFile *inputfile_qcd = TFile::Open("prospects/allQCD.root");
	TFile *inputfile_vv = TFile::Open("prospects/allVV.root");
	TFile *inputfile_wjets = TFile::Open("prospects/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");

	//loading needed variables
	Double_t norm = 1.;
	Double_t maximum = 0.;
	Double_t minimum = 10000.;
	Double_t temp_maximum = 0.;
	Double_t temp_minimum = 10000.;

	//woring on plots
	TH1F* h_dijetinvariantmass_sig = ((TH1F*)(inputfile_sig->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_sig->Scale(norm/h_dijetinvariantmass_sig->Integral("width"));
	h_dijetinvariantmass_sig->SetTitle("");
	h_dijetinvariantmass_sig->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_sig->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_sig->SetFillColor(7);
	h_dijetinvariantmass_sig->SetFillStyle(3001);
	h_dijetinvariantmass_sig->SetStats(kFALSE);
	temp_maximum = h_dijetinvariantmass_sig->GetMaximum();
	temp_minimum = h_dijetinvariantmass_sig->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_qcd = ((TH1F*)(inputfile_qcd->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_qcd->SetLineColor(kRed);
	//h_dijetinvariantmass_qcd->SetLineWidth(1.5);
	h_dijetinvariantmass_qcd->Scale(norm/h_dijetinvariantmass_qcd->Integral("width"));
	h_dijetinvariantmass_qcd->SetTitle("");
	h_dijetinvariantmass_qcd->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_qcd->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_qcd->SetStats(kFALSE);
	temp_maximum = h_dijetinvariantmass_qcd->GetMaximum();
	temp_minimum = h_dijetinvariantmass_qcd->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	TH1F* h_dijetinvariantmass_vv = ((TH1F*)(inputfile_vv->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_vv->SetLineColor(kBlack);
	h_dijetinvariantmass_vv->Scale(norm/h_dijetinvariantmass_vv->Integral("width"));
	h_dijetinvariantmass_vv->SetTitle("");
	h_dijetinvariantmass_vv->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_vv->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_vv->SetStats(kFALSE);
	temp_maximum = h_dijetinvariantmass_vv->GetMaximum();
	temp_minimum = h_dijetinvariantmass_vv->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	
	TH1F* h_dijetinvariantmass_wjets = ((TH1F*)(inputfile_wjets->Get(("demo/baselineObjectSelection/" + variable).c_str())));
	h_dijetinvariantmass_wjets->SetLineColor(kBlue);
	h_dijetinvariantmass_wjets->Scale(norm/h_dijetinvariantmass_wjets->Integral("width"));
	h_dijetinvariantmass_wjets->SetTitle("");
	h_dijetinvariantmass_wjets->GetXaxis()->SetTitle((label).c_str());
	h_dijetinvariantmass_wjets->GetYaxis()->SetTitle("a.u.");
	h_dijetinvariantmass_wjets->SetStats(kFALSE);
	temp_maximum = h_dijetinvariantmass_wjets->GetMaximum();
	temp_minimum = h_dijetinvariantmass_wjets->GetMinimum();
	if (temp_maximum > maximum) maximum = temp_maximum;
	if ((temp_minimum < minimum) && (temp_minimum != 0.) && ((temp_minimum > 0.00001)) ) minimum = temp_minimum;

	//Setting proper Y axis range
	maximum = maximum + 0.3 * maximum;
	h_dijetinvariantmass_sig->GetYaxis()->SetRangeUser(minimum,maximum);		
	h_dijetinvariantmass_qcd->GetYaxis()->SetRangeUser(minimum,maximum);		
	h_dijetinvariantmass_vv->GetYaxis()->SetRangeUser(minimum,maximum);		
	h_dijetinvariantmass_wjets->GetYaxis()->SetRangeUser(minimum,maximum);		

	//defining legend	
	TLegend* leg = new TLegend(0.58,0.7,0.9,0.9);
	leg->AddEntry(h_dijetinvariantmass_sig,"#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}jj","f");
	leg->AddEntry(h_dijetinvariantmass_qcd,"QCD","l");
	leg->AddEntry(h_dijetinvariantmass_vv,"VVjj","l");
	leg->AddEntry(h_dijetinvariantmass_vv,"W+jets","l");

	//plotting on canvas and writing to disk
	TCanvas *my_canvas = new TCanvas("my_canvas","my_canvas",500,500);
	my_canvas->cd();

	gPad->SetLogy();
	h_dijetinvariantmass_sig->Draw("HIST");
	h_dijetinvariantmass_qcd->Draw("HIST same");
	h_dijetinvariantmass_vv->Draw("HIST same");
	h_dijetinvariantmass_wjets->Draw("HIST same");
	leg->Draw();

	my_canvas->Print((variable + "_prospects13tev.pdf").c_str());
	my_canvas->Close();

}

//main function
void prospects_13Tev (){
	var_plot("h_dijetinvariantmass", "M_{(jet,jet)} [GeV]");
	var_plot("h_tau1pt", "P_{t}(#tau_{1}) [GeV]");
	var_plot("h_tau2pt", "P_{t}(#tau_{2}) [GeV]");
	var_plot("h_jet1pt", "P_{t}(jet_{1}) [GeV]");
	var_plot("h_met", "#slash{E}_{T} [GeV]");
}


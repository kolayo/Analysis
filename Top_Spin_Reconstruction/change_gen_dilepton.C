/*
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/orcun/Delphes-3.4.2
root
.include /home/orcun/Delphes-3.4.2
.include /home/orcun/Delphes-3.4.2/external
*/

#include "TH1.h"
#include "TSystem.h"
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
//#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#endif

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TBranch.h>
#include <TMath.h>
#include <TTree.h>

using namespace std;



void change_root_dilepton()
{

//	input=TFile::Open("evt_out_0-11_ntuple.root","read");
//	input_tree = static_cast<TTree*>(input->Get("delphTree"));
//	ExRootTreeReader *treeReader = new ExRootTreeReader("delphTree");

	gSystem->Load("libDelphes");
//	TChain chain("delphTree");
	TChain *chain = new TChain("Delphes");
        chain->Add("events_PYTHIA8_0_1.root");
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t nEntries = treeReader->GetEntries();

	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchMET = treeReader->UseBranch("MissingET");
	TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");

TClonesArray *branchGen = treeReader->UseBranch("Particle");
TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
TClonesArray *branchGenMET = treeReader->UseBranch("GenMissingET");
        TFile *output = TFile::Open("reduced_events/Gen-Level/leppt10-3_jetpt10-5/ttbar-gen_cut-dilepton-output_1.root", "RECREATE");
	TTree *out_tree = new TTree("nominal","nominal");


	std::vector<float> gen_lepton_pt;
	std::vector<float> gen_lepton_eta;
	std::vector<float> gen_lepton_phi;
	std::vector<float> gen_lepton_e;
	std::vector<float> gen_lepton_charge;

	std::vector<float> lepton_pt;
	std::vector<float> lepton_eta;
	std::vector<float> lepton_phi;
	std::vector<float> lepton_e;
	std::vector<float> lepton_charge;
	float met_met=0;
	float met_phi=0;
	float sumet=0;
	std::vector<char> lepton_is_e_mu;
	std::vector<float> jet_pt;
	std::vector<float> jet_eta;
	std::vector<float> jet_phi;
	std::vector<float> jet_e;
	std::vector<float> jet_btag_weight;
	std::vector<char> jet_has_btag;

  std::vector<double> klf_costheta_k1;
  std::vector<double> klf_costheta_k2;
  std::vector<double> klf_costheta_r1;
  std::vector<double> klf_costheta_r2;
  std::vector<double> klf_costheta_n1;
  std::vector<double> klf_costheta_n2;
  std::vector<double> klf_costheta_k1_star;
  std::vector<double> klf_costheta_k2_star;
  std::vector<double> klf_costheta_r1_star;
  std::vector<double> klf_costheta_r2_star;
  std::vector<double> klf_C_kk;
  std::vector<double> klf_C_rr;
  std::vector<double> klf_C_nn;
  std::vector<double> klf_C_rk_plus_kr;
  std::vector<double> klf_C_rk_minus_kr;
  std::vector<double> klf_C_nr_plus_rn;
  std::vector<double> klf_C_nr_minus_rn;
  std::vector<double> klf_C_nk_plus_kn;
  std::vector<double> klf_C_nk_minus_kn;
  std::vector<double> klf_cosphi_rest_ttbar;
  std::vector<double> klf_cosphi_lab;
  std::vector<double> klf_delta_phi;


	out_tree->Branch("gen_lepton_pt", &gen_lepton_pt);
	out_tree->Branch("gen_lepton_eta", &gen_lepton_eta);
	out_tree->Branch("gen_lepton_phi", &gen_lepton_phi);
	out_tree->Branch("gen_lepton_e", &gen_lepton_e);
	out_tree->Branch("gen_lepton_charge", &gen_lepton_charge);
	out_tree->Branch("lepton_pt", &lepton_pt);
	out_tree->Branch("lepton_eta", &lepton_eta);
	out_tree->Branch("lepton_phi", &lepton_phi);
	out_tree->Branch("lepton_e", &lepton_e);
	out_tree->Branch("lepton_charge", &lepton_charge);
	out_tree->Branch("met_met", &met_met);
	out_tree->Branch("met_phi", &met_phi);
	out_tree->Branch("sumet", &sumet);
	out_tree->Branch("lepton_is_e_mu", &lepton_is_e_mu);
	out_tree->Branch("jet_pt", &jet_pt);
	out_tree->Branch("jet_eta", &jet_eta);
	out_tree->Branch("jet_phi", &jet_phi);
	out_tree->Branch("jet_e", &jet_e);
	out_tree->Branch("jet_btag_weight", &jet_btag_weight);
	out_tree->Branch("jet_has_btag", &jet_has_btag);

  out_tree->Branch("klf_costheta_k1", &klf_costheta_k1);
  out_tree->Branch("klf_costheta_r1", &klf_costheta_r1);
  out_tree->Branch("klf_costheta_n1", &klf_costheta_n1);
  out_tree->Branch("klf_costheta_k2", &klf_costheta_k2);
  out_tree->Branch("klf_costheta_r2", &klf_costheta_r2);
  out_tree->Branch("klf_costheta_n2", &klf_costheta_n2);
  out_tree->Branch("klf_costheta_k1_star", &klf_costheta_k1_star);
  out_tree->Branch("klf_costheta_r1_star", &klf_costheta_r1_star);
  out_tree->Branch("klf_costheta_k2_star", &klf_costheta_k2_star);
  out_tree->Branch("klf_costheta_r2_star", &klf_costheta_r2_star);
  out_tree->Branch("klf_C_kk", &klf_C_kk);
  out_tree->Branch("klf_C_rr", &klf_C_rr);
  out_tree->Branch("klf_C_nn", &klf_C_nn);
  out_tree->Branch("klf_C_rk_plus_kr", &klf_C_rk_plus_kr);
  out_tree->Branch("klf_C_rk_minus_kr", &klf_C_rk_minus_kr);
  out_tree->Branch("klf_C_nr_plus_rn", &klf_C_nr_plus_rn);
  out_tree->Branch("klf_C_nr_minus_rn", &klf_C_nr_minus_rn);
  out_tree->Branch("klf_C_nk_plus_kn", &klf_C_nk_plus_kn);
  out_tree->Branch("klf_C_nk_minus_kn", &klf_C_nk_minus_kn);
  out_tree->Branch("klf_cosphi_rest_ttbar", &klf_cosphi_rest_ttbar);
  out_tree->Branch("klf_cosphi_lab", &klf_cosphi_lab);
  out_tree->Branch("klf_delta_phi", &klf_delta_phi);
// ******************************************************************************************************
  TH1F* hjet_pt     = new TH1F("hjet_pt","hjet_pt",50,0,500);
  TH1F* hjet_eta     = new TH1F("hjet_eta","hjet_eta",50,-4,4);
  TH1F* hjet_phi     = new TH1F("hjet_phi","hjet_phi",50,-4,4);
  TH1F* hMET     = new TH1F("hMET","hMET",50,0,500);
  TH1F* hb_jet_pt     = new TH1F("hb_jet_pt","hb_jet_pt",50,0,500);
  TH1F* hb_jet_eta     = new TH1F("hb_jet_eta","hb_jet_eta",50,-4,4);
  TH1F* hb_jet_phi     = new TH1F("hb_jet_phi","hb_jet_phi",50,-4,4);
  TH1F* hlep_pt     = new TH1F("hlep_pt","hlep_pt",50,0,500);
  TH1F* hlep_eta     = new TH1F("hlep_eta","hlep_eta",50,-4,4);
  TH1F* hlep_phi     = new TH1F("hlep_phi","hlep_phi",50,-4,4);
  TH1F* hcostheta_k1     = new TH1F("hcostheta_k1","B^{k}_{1}",6,-1,1);
  TH1F* hcostheta_r1     = new TH1F("hcostheta_r1","B^{r}_{1}",6,-1,1);
  TH1F* hcostheta_n1     = new TH1F("hcostheta_n1","B^{n}_{1}",6,-1,1);
  TH1F* hcostheta_k2     = new TH1F("hcostheta_k2","B^{k}_{2}",6,-1,1);
  TH1F* hcostheta_r2     = new TH1F("hcostheta_r2","B^{r}_{2}",6,-1,1);
  TH1F* hcostheta_n2     = new TH1F("hcostheta_n2","B^{n}_{2}",6,-1,1);
  TH1F* hcostheta_k1_star     = new TH1F("hcostheta_k1_star","B^{k*}_{1}",6,-1,1);
  TH1F* hcostheta_r1_star     = new TH1F("hcostheta_r1_star","B^{r*}_{1}",6,-1,1);
  TH1F* hcostheta_k2_star     = new TH1F("hcostheta_k2_star","B^{k*}_{2}",6,-1,1);
  TH1F* hcostheta_r2_star     = new TH1F("hcostheta_r2_star","B^{r*}_{2}",6,-1,1);
  TH1F* hC_kk     = new TH1F("hC_kk","C_{kk}",6,-1,1);
  TH1F* hC_rr     = new TH1F("hC_rr","C_{rr}",6,-1,1);
  TH1F* hC_nn     = new TH1F("hC_nn","C_{nn}",6,-1,1);
  TH1F* hC_rk_plus_kr     = new TH1F("hC_rk_plus_kr","C_{rk}+C_{kr}",6,-1,1);
  TH1F* hC_rk_minus_kr     = new TH1F("hC_rk_minus_kr","C_{rk}-C_{kr}",6,-1,1);
  TH1F* hC_nr_plus_rn     = new TH1F("hC_nr_plus_rn","C_{nr}+C_{rn}",6,-1,1);
  TH1F* hC_nr_minus_rn     = new TH1F("hC_nr_minus_rn","C_{nr}-C_{rn}",6,-1,1);
  TH1F* hC_nk_plus_kn     = new TH1F("hC_nk_plus_kn","C_{nk}+C_{kn}",6,-1,1);
  TH1F* hC_nk_minus_kn     = new TH1F("hC_nk_minus_kn","C_{nk}-C_{kn}",6,-1,1);
  TH1F* hcosphi_rest_ttbar     = new TH1F("hcosphi_rest_ttbar","cos#varphi",6,-1,1);
  TH1F* hcosphi_lab     = new TH1F("hcosphi_lab","cos#varphi_{lab}",6,-1,1);
  TH1F* hdelta_phi     = new TH1F("hdelta_phi","|#Delta#phi_{ll}|",6,0,3.14159);

// ******************************************************************************************************


int lep_say = 0;
int say_top = 0;
int say_evt = 0;
int real_event = 0;
	int Nel = 0;
	int Nmu = 0;
	int numberofjets=0;
int Nparticle=0;
int say=0;
	std::cout << "Started looping over " << nEntries << " entries" << std::endl;
	Long64_t ievent=0;
	for (ievent = 0; ievent < nEntries; ++ievent) {
		gen_lepton_pt.clear();
		gen_lepton_eta.clear();
		gen_lepton_phi.clear();
		gen_lepton_e.clear();
		gen_lepton_charge.clear();
		lepton_pt.clear();
		lepton_eta.clear();
		lepton_phi.clear();
		lepton_e.clear();
		lepton_charge.clear();
		lepton_is_e_mu.clear();
		jet_pt.clear();
		jet_eta.clear();
		jet_phi.clear();
		jet_e.clear();
		jet_btag_weight.clear();
		jet_has_btag.clear();
		met_met=met_phi=sumet=0;

    klf_costheta_k1.clear();
    klf_costheta_r1.clear();
    klf_costheta_n1.clear();
    klf_costheta_k2.clear();
    klf_costheta_r2.clear();
    klf_costheta_n2.clear();
    klf_costheta_k1_star.clear();
    klf_costheta_r1_star.clear();
    klf_costheta_k2_star.clear();
    klf_costheta_r2_star.clear();
    klf_C_kk.clear();
    klf_C_rr.clear();
    klf_C_nn.clear();
    klf_C_rk_plus_kr.clear();
    klf_C_rk_minus_kr.clear();
    klf_C_nr_plus_rn.clear();
    klf_C_nr_minus_rn.clear();
    klf_C_nk_plus_kn.clear();
    klf_C_nk_minus_kn.clear();
    klf_cosphi_rest_ttbar.clear();
    klf_cosphi_lab.clear();
    klf_delta_phi.clear();




		treeReader->ReadEntry(ievent);
		Nel = branchElectron->GetEntries();
		Nmu = branchMuon->GetEntries();
		numberofjets = branchJet->GetEntries();

Nparticle = branchGen->GetEntries();
//cout << "number of particle :  " << Nparticle << endl;
cout << "number of event: " << ievent << endl;
TLorentzVector Top_lab;
TLorentzVector Topbar_lab;
TLorentzVector Positive_lep_lab;
TLorentzVector Negative_lep_lab;
GenParticle *top;
GenParticle *topbar;
GenParticle *lep_pos;
GenParticle *lep_neg;
double Bottom_Pt = -9999;
double Bottom_Eta = -9999;
double Bottom_Phi = -9999;
double Bottombar_Pt = -9999;
double Bottombar_Eta = -9999;
double Bottombar_Phi = -9999;
double Positive_lep_Pt = -9999;
double Positive_lep_Eta = -9999;
double Negative_lep_Pt = -9999;
double Negative_lep_Eta = -9999;
double Positive_lep_Phi = -9999;
double Negative_lep_Phi = -9999;
int pos_lep = 0;
int neg_lep = 0;
int dilepton = 0;
int top_index = 0;
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//Selection of Correct Leptons and Top Quarks//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-/
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
for (int i=0; i<Nparticle; ++i){
	GenParticle *genlep = (GenParticle*) branchGen->At(i);
//cout << "(INDEX PID D1 D2 STATUS MASS PX PY PZ E): " << i << " " << genlep->PID << " " << genlep->D1 << " " << genlep->D2 << " " << genlep->Status << " " << genlep->Mass << " " << genlep->Px << " " << genlep->Py << " " << genlep->Pz <<" " << genlep->E << endl;
	if (genlep->M1 == -1) continue;
	GenParticle *mother_genlep = (GenParticle*) branchGen->At(genlep->M1);
	if (mother_genlep->M1 == -1) continue;
	GenParticle *gmother_genlep = (GenParticle*) branchGen->At(mother_genlep->M1);
	top_index = mother_genlep->M1;
//	if (!(abs(mother_genlep->PID) == 24) || !(abs(genlep->PID) == 11 || abs(genlep->PID) == 13 || abs(genlep->PID) == 15)) continue;
	if (!(abs(mother_genlep->PID) == 24) || !(abs(genlep->PID) == 11 || abs(genlep->PID) == 13)) continue;
	if (genlep->Status != 22 && genlep->Status != 23 && genlep->Status != 1) continue;   //****************************************************************************

	while (abs(gmother_genlep->PID) == 24) {
		top_index = gmother_genlep->M1;
		gmother_genlep = (GenParticle*) branchGen->At(gmother_genlep->M1);
//cout << "GMOTHER_GENLEP: " << " PID: " << gmother_genlep->PID << " D1: " << gmother_genlep->D1 << " D2: " << gmother_genlep->D2 << " STATUS: " << gmother_genlep->Status << " MASS: " << gmother_genlep->Mass << endl;
//		GenParticle *dau1_gmother_genlep = (GenParticle*) branchGen->At(gmother_genlep->D1);
//		GenParticle *dau2_gmother_genlep = (GenParticle*) branchGen->At(gmother_genlep->D2);
//cout << "DAU1_GMOTHER_GENLEP: " << " PID: " << dau1_gmother_genlep->PID << " D1: " << dau1_gmother_genlep->D1 << " D2: " << dau1_gmother_genlep->D2 << " STATUS: " << dau1_gmother_genlep->Status << " MASS: " << dau1_gmother_genlep->Mass << endl;
//cout << "DAU2_GMOTHER_GENLEP: " << " PID: " << dau2_gmother_genlep->PID << " D1: " << dau2_gmother_genlep->D1 << " D2: " << dau2_gmother_genlep->D2 << " STATUS: " << dau2_gmother_genlep->Status << " MASS: " << dau2_gmother_genlep->Mass << endl;
	}
        if (gmother_genlep->Status != 62 ) continue;
	if (gmother_genlep->PID == 6) {
		top = (GenParticle*) branchGen->At(top_index);
//cout << "TOP INDEX AND PID: " << top_index << " " << top->PID << endl;
	}
	if (gmother_genlep->PID == -6) {
		topbar = (GenParticle*) branchGen->At(top_index);
//cout << "TOPBAR INDEX AND PID: " << top_index << " " << topbar->PID << endl;
	}
	
	if (genlep->PID < 0) {
		lep_pos = (GenParticle*) branchGen->At(i);
	}
	if (genlep->PID > 0) {
		lep_neg = (GenParticle*) branchGen->At(i);
	}
//cout << "LEPTON MOTHER GMOTHER INDEX AND PID: " << i << " " << genlep->M1 << " " << mother_genlep->M1 << " " << genlep->PID << " " << mother_genlep->PID << " " << gmother_genlep->PID << endl;


dilepton++;
}
if (dilepton != 2) continue;
GenParticle *bottom1 = (GenParticle*) branchGen->At(top->D1);                            // D1 is bottom quark in MadGraph, but D2 is bottom in Sherpa.
GenParticle *bottom2 = (GenParticle*) branchGen->At(topbar->D1);                         // D1 is bottom quark in MadGraph, but D2 is bottom in Sherpa.
cout << bottom1->PID << " " << bottom2->PID << endl;
if (lep_pos->PT <= 10 || abs(lep_pos->Eta) > 3 || bottom1->PT <= 10 || abs(bottom1->Eta) > 5) continue;
if (lep_neg->PT <= 10 || abs(lep_neg->Eta) > 3 || bottom2->PT <= 10 || abs(bottom2->Eta) > 5) continue;
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//Selection of Correct Leptons and Top Quarks//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-/
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
Top_lab.SetPxPyPzE(top->Px, top->Py, top->Pz, top->E);
Topbar_lab.SetPxPyPzE(topbar->Px, topbar->Py, topbar->Pz, topbar->E);
Positive_lep_lab.SetPxPyPzE(lep_pos->Px, lep_pos->Py, lep_pos->Pz, lep_pos->E);
Negative_lep_lab.SetPxPyPzE(lep_neg->Px, lep_neg->Py, lep_neg->Pz, lep_neg->E);
TVector3 Positive_lep_lab_v3(Positive_lep_lab.Px(), Positive_lep_lab.Py(), Positive_lep_lab.Pz());
TVector3 Negative_lep_lab_v3(Negative_lep_lab.Px(), Negative_lep_lab.Py(), Negative_lep_lab.Pz());
//************************************************************************************************
Positive_lep_Pt = lep_pos->PT;
Positive_lep_Eta = lep_pos->Eta;
Positive_lep_Phi = lep_pos->Phi;
Negative_lep_Pt = lep_neg->PT;
Negative_lep_Eta = lep_neg->Eta;
Negative_lep_Phi = lep_neg->Phi;
hb_jet_pt->Fill(bottom1->PT);
hb_jet_eta->Fill(bottom1->Eta);
hb_jet_phi->Fill(bottom1->Phi);
hb_jet_pt->Fill(bottom2->PT);
hb_jet_eta->Fill(bottom2->Eta);
hb_jet_phi->Fill(bottom2->Phi);
hlep_pt->Fill(lep_pos->PT);
hlep_eta->Fill(lep_pos->Eta);
hlep_phi->Fill(lep_pos->Phi);
hlep_pt->Fill(lep_neg->PT);
hlep_eta->Fill(lep_neg->Eta);
hlep_phi->Fill(lep_neg->Phi);
int Ngen_jets = 0;
Ngen_jets = branchGenJet->GetEntries();
for (int j=0; j<Ngen_jets; ++j) {
	Jet *genjet = (Jet*) branchGenJet->At(j);
	hjet_pt->Fill(genjet->PT);
	hjet_eta->Fill(genjet->Eta);
	hjet_phi->Fill(genjet->Phi);
}
MissingET *genmet = (MissingET*) branchGenMET->At(0);
hMET->Fill(genmet->MET);
//************************************************************************************************
real_event++;
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//Reconstruct Spin Observables//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
TVector3 TopTopbar_lab(Top_lab.Px()+Topbar_lab.Px(), Top_lab.Py()+Topbar_lab.Py(), Top_lab.Pz()+Topbar_lab.Pz());
TVector3 beta1(-TopTopbar_lab.Px()/(Top_lab.E()+Topbar_lab.E()), -TopTopbar_lab.Py()/(Top_lab.E()+Topbar_lab.E()), -TopTopbar_lab.Pz()/(Top_lab.E()+Topbar_lab.E()));
TLorentzVector Top_ttbar_zmf;
Top_ttbar_zmf.SetPxPyPzE(Top_lab.Px(), Top_lab.Py(), Top_lab.Pz(), Top_lab.E());
Top_ttbar_zmf.Boost(beta1);
TVector3 Top_ttbar_zmf_v3(Top_ttbar_zmf.Px(), Top_ttbar_zmf.Py(), Top_ttbar_zmf.Pz());
TLorentzVector Topbar_ttbar_zmf;
Topbar_ttbar_zmf.SetPxPyPzE(Topbar_lab.Px(), Topbar_lab.Py(), Topbar_lab.Pz(), Topbar_lab.E());
Topbar_ttbar_zmf.Boost(beta1);
TVector3 Topbar_ttbar_zmf_v3(Topbar_ttbar_zmf.Px(), Topbar_ttbar_zmf.Py(), Topbar_ttbar_zmf.Pz());
TVector3 p_beam(0.0, 0.0, 1.0);
TVector3 top_Pk;
TVector3 top_Pr;
TVector3 top_Pn;
TVector3 topbar_Pk;
TVector3 topbar_Pr;
TVector3 topbar_Pn;
top_Pk.SetXYZ(Top_ttbar_zmf_v3.Unit().Px(), Top_ttbar_zmf_v3.Unit().Py(), Top_ttbar_zmf_v3.Unit().Pz());
topbar_Pk = -1*top_Pk;
double y = p_beam.Dot(top_Pk);                                                            // cos(theta) between k and p
double r = sqrt(1-y*y);
if (y>=0) {
	top_Pr = 1/r*(p_beam-y*top_Pk);
	top_Pn = 1/r*(p_beam.Cross(top_Pk));
	topbar_Pr = -1*top_Pr;
	topbar_Pn = -1*top_Pn;
}
if (y<0) {
	top_Pr = -1/r*(p_beam-y*top_Pk);
	top_Pn = -1/r*(p_beam.Cross(top_Pk));
	topbar_Pr = -1*top_Pr;
	topbar_Pn = -1*top_Pn;
}
TLorentzVector lep1_top_cmf;                                                                            // boosted (parent CM frame) lep1 (positively charged lepton)
TLorentzVector lep2_topbar_cmf;                                                                         // boosted (parent CM frame) lep2 (negatively charged lepton)
lep1_top_cmf.SetPxPyPzE(Positive_lep_lab.Px(), Positive_lep_lab.Py(), Positive_lep_lab.Pz(), Positive_lep_lab.E());
lep2_topbar_cmf.SetPxPyPzE(Negative_lep_lab.Px(), Negative_lep_lab.Py(), Negative_lep_lab.Pz(), Negative_lep_lab.E());
lep1_top_cmf.Boost(beta1);
lep2_topbar_cmf.Boost(beta1);
TVector3 beta_top(-Top_ttbar_zmf.Px()/Top_ttbar_zmf.E(), -Top_ttbar_zmf.Py()/Top_ttbar_zmf.E(), -Top_ttbar_zmf.Pz()/Top_ttbar_zmf.E());
TVector3 beta_topbar(-Topbar_ttbar_zmf.Px()/Topbar_ttbar_zmf.E(), -Topbar_ttbar_zmf.Py()/Topbar_ttbar_zmf.E(), -Topbar_ttbar_zmf.Pz()/Topbar_ttbar_zmf.E());
lep1_top_cmf.Boost(beta_top);
lep2_topbar_cmf.Boost(beta_topbar);
TVector3 lep1_top_cmf_v3(lep1_top_cmf.Px(), lep1_top_cmf.Py(), lep1_top_cmf.Pz());
TVector3 lep2_topbar_cmf_v3(lep2_topbar_cmf.Px(), lep2_topbar_cmf.Py(), lep2_topbar_cmf.Pz());

          double costheta_k1 = -9999;
            double costheta_r1 = -9999;
            double costheta_n1 = -9999;
            double costheta_k2 = -9999;
            double costheta_r2 = -9999;
            double costheta_n2 = -9999;
            double costheta_k1_star = -9999;
            double costheta_r1_star = -9999;
            double costheta_k2_star = -9999;
            double costheta_r2_star = -9999;
            double C_kk = -9999;
            double C_rr = -9999;
            double C_nn = -9999;
            double C_rk_plus_kr = -9999;
            double C_rk_minus_kr = -9999;
            double C_nr_plus_rn = -9999;
            double C_nr_minus_rn = -9999;
            double C_nk_plus_kn = -9999;
            double C_nk_minus_kn = -9999;
            double cosphi_rest_ttbar = -9999;
            double cosphi_lab = -9999;
            double delta_phi = -9999;
      costheta_k1 = lep1_top_cmf_v3.Unit().Dot(top_Pk);         // angle between lep1 (positively charged lepton) and top k axis in ttbar zmf
      costheta_k2 = lep2_topbar_cmf_v3.Unit().Dot(topbar_Pk);   // angle between lep2 (negatively charged lepton) and topbar k axis in ttbar zmf
      costheta_r1 = lep1_top_cmf_v3.Unit().Dot(top_Pr);         // angle between lep1 (positively charged lepton) and top r axis in ttbar zmf
      costheta_r2 = lep2_topbar_cmf_v3.Unit().Dot(topbar_Pr);   // angle between lep2 (negatively charged lepton) and topbar r axis in ttbar zmf
      costheta_n1 = lep1_top_cmf_v3.Unit().Dot(top_Pn);         // angle between lep1 (positively charged lepton) and top n axis in ttbar zmf
      costheta_n2 = lep2_topbar_cmf_v3.Unit().Dot(topbar_Pn);   // angle between lep2 (negatively charged lepton) and topbar n axis in ttbar zmf
      double rapidity_top = atanh(Top_lab.Pz()/Top_lab.E());
      double rapidity_topbar = atanh(Topbar_lab.Pz()/Topbar_lab.E());
      double delta_rapidity = abs(rapidity_top)-abs(rapidity_topbar);
      if (delta_rapidity>=0) {
        costheta_k1_star = costheta_k1;
        costheta_r1_star = costheta_r1;
        costheta_k2_star = costheta_k2;
        costheta_r2_star = costheta_r2;
      }
      if (delta_rapidity<0) {
        costheta_k1_star = -costheta_k1;
        costheta_r1_star = -costheta_r1;
        costheta_k2_star = -costheta_k2;
        costheta_r2_star = -costheta_r2;
      }
      C_kk = costheta_k1*costheta_k2;
      C_rr = costheta_r1*costheta_r2;
      C_nn = costheta_n1*costheta_n2;
      C_rk_plus_kr = costheta_r1*costheta_k2+costheta_k1*costheta_r2;
      C_rk_minus_kr = costheta_r1*costheta_k2-costheta_k1*costheta_r2;
      C_nr_plus_rn = costheta_n1*costheta_r2+costheta_r1*costheta_n2;
      C_nr_minus_rn = costheta_n1*costheta_r2-costheta_r1*costheta_n2;
      C_nk_plus_kn = costheta_n1*costheta_k2+costheta_k1*costheta_n2;
      C_nk_minus_kn = costheta_n1*costheta_k2-costheta_k1*costheta_n2;
      cosphi_rest_ttbar = lep1_top_cmf_v3.Unit().Dot(lep2_topbar_cmf_v3.Unit());
      cosphi_lab = Positive_lep_lab_v3.Unit().Dot(Negative_lep_lab_v3.Unit());
      if (Positive_lep_Phi >= Negative_lep_Phi) {
        if ((Positive_lep_Phi-Negative_lep_Phi)>acos(-1)) {
          delta_phi = 2*acos(-1)-(Positive_lep_Phi-Negative_lep_Phi);
        }
        if ((Positive_lep_Phi-Negative_lep_Phi)<=acos(-1)) {
          delta_phi = (Positive_lep_Phi-Negative_lep_Phi);
        }
      }
      if (Positive_lep_Phi < Negative_lep_Phi) {
        if ((Negative_lep_Phi-Positive_lep_Phi)>acos(-1)) {
          delta_phi = 2*acos(-1)-(Negative_lep_Phi-Positive_lep_Phi);
        }
        if ((Negative_lep_Phi-Positive_lep_Phi)<=acos(-1)) {
          delta_phi = (Negative_lep_Phi-Positive_lep_Phi);
        }
      }
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//Reconstruct Spin Observables//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//

      klf_costheta_k1.emplace_back(costheta_k1);
      klf_costheta_r1.emplace_back(costheta_r1);
      klf_costheta_n1.emplace_back(costheta_n1);
      klf_costheta_k2.emplace_back(costheta_k2);
      klf_costheta_r2.emplace_back(costheta_r2);
      klf_costheta_n2.emplace_back(costheta_n2);
      klf_costheta_k1_star.emplace_back(costheta_k1_star);
      klf_costheta_r1_star.emplace_back(costheta_r1_star);
      klf_costheta_k2_star.emplace_back(costheta_k2_star);
      klf_costheta_r2_star.emplace_back(costheta_r2_star);
      klf_C_kk.emplace_back(C_kk);
      klf_C_rr.emplace_back(C_rr);
      klf_C_nn.emplace_back(C_nn);
      klf_C_rk_plus_kr.emplace_back(C_rk_plus_kr);
      klf_C_rk_minus_kr.emplace_back(C_rk_minus_kr);
      klf_C_nr_plus_rn.emplace_back(C_nr_plus_rn);
      klf_C_nr_minus_rn.emplace_back(C_nr_minus_rn);
      klf_C_nk_plus_kn.emplace_back(C_nk_plus_kn);
      klf_C_nk_minus_kn.emplace_back(C_nk_minus_kn);
      klf_cosphi_rest_ttbar.emplace_back(cosphi_rest_ttbar);
      klf_cosphi_lab.emplace_back(cosphi_lab);
      klf_delta_phi.emplace_back(delta_phi);

      hcostheta_k1->Fill(costheta_k1);
      hcostheta_r1->Fill(costheta_r1);
      hcostheta_n1->Fill(costheta_n1);
      hcostheta_k2->Fill(costheta_k2);
      hcostheta_r2->Fill(costheta_r2);
      hcostheta_n2->Fill(costheta_n2);
      hcostheta_k1_star->Fill(costheta_k1_star);
      hcostheta_r1_star->Fill(costheta_r1_star);
      hcostheta_k2_star->Fill(costheta_k2_star);
      hcostheta_r2_star->Fill(costheta_r2_star);
      hC_kk->Fill(C_kk);
      hC_rr->Fill(C_rr);
      hC_nn->Fill(C_nn);
      hC_rk_plus_kr->Fill(C_rk_plus_kr);
      hC_rk_minus_kr->Fill(C_rk_minus_kr);
      hC_nr_plus_rn->Fill(C_nr_plus_rn);
      hC_nr_minus_rn->Fill(C_nr_minus_rn);
      hC_nk_plus_kn->Fill(C_nk_plus_kn);
      hC_nk_minus_kn->Fill(C_nk_minus_kn);
      hcosphi_rest_ttbar->Fill(cosphi_rest_ttbar);
      hcosphi_lab->Fill(cosphi_lab);
      hdelta_phi->Fill(delta_phi);

	out_tree->Fill();
	}
//	std::cout << "semileptonic events " << say  << std::endl;
cout << "DEBUG TOTAL : " << say_evt << endl;
cout << "total top/anti-top : " << say_top << endl;
cout << "total event: " << real_event << endl;
cout << "total ilk lepton: " << lep_say << endl;
output->cd();
out_tree->Write();
//out_tree->Print();
  hjet_pt->Write();
  hjet_eta->Write();
  hjet_phi->Write();
  hMET->Write();
  hb_jet_pt->Write();
  hb_jet_eta->Write();
  hb_jet_phi->Write();
  hlep_pt->Write();
  hlep_eta->Write();
  hlep_phi->Write();
  hcostheta_k1->Write();
  hcostheta_r1->Write();
  hcostheta_n1->Write();
  hcostheta_k2->Write();
  hcostheta_r2->Write();
  hcostheta_n2->Write();
  hcostheta_k1_star->Write();
  hcostheta_r1_star->Write();
  hcostheta_k2_star->Write();
  hcostheta_r2_star->Write();
  hC_kk->Write();
  hC_rr->Write();
  hC_nn->Write();
  hC_rk_plus_kr->Write();
  hC_rk_minus_kr->Write();
  hC_nr_plus_rn->Write();
  hC_nr_minus_rn->Write();
  hC_nk_plus_kn->Write();
  hC_nk_minus_kn->Write();
  hcosphi_rest_ttbar->Write();
  hcosphi_lab->Write();
  hdelta_phi->Write();


}

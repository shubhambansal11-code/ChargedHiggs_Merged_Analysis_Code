#ifndef EventLoop_h
#define EventLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "string"
#include "TLorentzVector.h"
#include "vector"
#include <iostream>

#include "utilis/Chi2_minimization.h"	
#include "utilis/NeutrinoBuilder.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#include "TH1Fs/TH1Fs.h"
#include "TMVA/Reader.h"

class EventLoop {
public :
   TTree           *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   // Declaration of leaf types
   Int_t           nJet;
   Int_t           mcChannelNumber;
   Int_t           TopHeavyFlavorFilterFlag;
   Int_t           nMuons;
   Int_t           nElectrons;
   Long_t          eventNumber;
   Int_t           nPrimaryVtx;
   Int_t           RunNumber;
   Int_t           randomRunNumber;
   Int_t           isTriggered;
   Int_t           nTaus;
   Int_t           HF_SimpleClassification;
   Int_t           jet_truthflav; 
   Float_t         Mu;
   Float_t         EventWeight;

   
   Float_t         weight_pileup;
   Float_t         weight_mc;
   Float_t         weight_leptonSF;
   Float_t         weight_trackjet_bTagSF_DL1r_Continuous;
   Float_t         weight_jvt;
   
   Float_t         met;
   Float_t         met_phi;
   std::vector<float>   *el_pt;
   std::vector<float>   *el_eta;
   std::vector<float>   *el_phi;
   std::vector<float>   *el_e;
   std::vector<float>   *el_charge;
   std::vector<float>   *mu_pt;
   std::vector<float>   *mu_eta;
   std::vector<float>   *mu_phi;
   std::vector<float>   *mu_e;
   std::vector<float>   *mu_charge;
   
   std::vector<float>   *ljet_m;
   std::vector<float>   *ljet_pt;
   std::vector<float>   *ljet_eta;
   std::vector<float>   *ljet_phi;
   std::vector<int>   *ljet_truthLabel;

   std::vector<float>   *ljet_D2;
   std::vector<float>   *ljet_C2;
   std::vector<float>   *ljet_Xbb2020v3_Higgs;
   std::vector<float>   *ljet_Xbb2020v3_QCD;
   std::vector<float>   *ljet_Xbb2020v3_Top;
   
   std::vector<float>   *tjet_e;
   std::vector<float>   *tjet_pt;
   std::vector<float>   *tjet_eta;
   std::vector<float>   *tjet_phi;
   std::vector<int>   *tjet_tagWeightBin_DL1r_Continuous;
  
   std::vector<float>   *signal_Jet_E;
   std::vector<float>   *signal_Jet_PT;
   std::vector<float>   *signal_Jet_Eta;
   std::vector<float>   *signal_Jet_Phi;
   std::vector<int>   *signal_Jet_truthflav;  
   std::vector<int>     *signal_Jet_tagWeightBin_DL1r_Continuous;

   std::vector<int>     *truth_pdgid;
   std::vector<int>     *truth_status;
   std::vector<float>   *truth_pt; 
   std::vector<float>   *truth_eta; 
   std::vector<float>   *truth_phi; 
   std::vector<float>   *truth_m;
   std::vector<long>   *truth_tthbb_info; 

   std::vector<TLorentzVector*> NeutrinoVec;
   map<TString, bool> pass_sel;

   // List of branches
   TBranch        *b_nJet;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_randomRunNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nPrimaryVtx;   //!
   TBranch        *b_nTaus;   //!
   TBranch        *b_HF_SimpleClassification; //!
   TBranch        *b_mcChannelNumber;
   TBranch        *b_TopHeavyFlavorFilterFlag;
   
   
   TBranch        *b_weight_pileup;   //!
   TBranch        *b_weight_mc;   //!
   TBranch        *b_weight_leptonSF;   //!
   TBranch        *b_weight_trackjet_bTagSF_DL1r_Continuous;   //!
   TBranch        *b_weight_jvt;   //! 
   

   TBranch        *b_MET;   //!
   TBranch        *b_MET_phi;   //!
  
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_signal_Jet_E;   //!
   TBranch        *b_signal_Jet_PT;   //!
   TBranch        *b_signal_Jet_Eta;   //!
   TBranch        *b_signal_Jet_Phi;   //!
   TBranch        *b_signal_Jet_truthflav;   //!   //only for MC
   TBranch        *b_signal_Jet_tagWeightBin_DL1r_Continuous;   //!

   TBranch        *b_ljet_m;   //!
   TBranch        *b_ljet_pt;   //!
   TBranch        *b_ljet_eta;   //!
   TBranch        *b_ljet_phi;   //!
   TBranch        *b_ljet_truthLabel;   //!   

   TBranch   *b_ljet_D2;  //!
   TBranch   *b_ljet_C2;  //!
   TBranch   *b_ljet_Xbb2020v3_Higgs; //!
   TBranch   *b_ljet_Xbb2020v3_QCD;  //!
   TBranch   *b_ljet_Xbb2020v3_Top;  //!
   
   TBranch   *b_tjet_e; //!
   TBranch   *b_tjet_pt; //!
   TBranch   *b_tjet_eta; //!
   TBranch   *b_tjet_phi; //!
   TBranch   *b_tjet_tagWeightBin_DL1r_Continuous; //!
   
   TBranch        *b_truth_pt;   //!
   TBranch        *b_truth_eta;   //!
   TBranch        *b_truth_phi;   //!
   TBranch        *b_truth_m;   //!
   TBranch        *b_truth_pdgid;   //! 
   TBranch        *b_truth_status;   //!
   TBranch        *b_truth_tthbb_info; //!  
   

   EventLoop(TTree *tree=0, TString sampleName="", TString MCDataPeriode="", TString ExpUncertaintyName="Nominal", TString WP="");
   void Write(TDirectory *dir, std::string dirname);
   void WriteTreesToFile(TFile *outFile);

   void	Sort_Jets(std::vector<TLorentzVector> *Jets, std::vector<int> *is_tagged);
   void SetTruthParticles();

   void SetJetVectors();
   void SetLeptonVectors();
   
   void dRjj();
   bool Findhad_tt();
   void Calc_dR_jl(float &dR_jl_min, float &dR_jl_max);
   double GetMwt();
   double GetTruthMass();
   
   bool FindFJetPair();
   bool PassEventSelectionBoosted();
   int GetBTagCategoryShort(int NTags_InHiggsJet, int NTags_OutsideHiggsJet);

   TLorentzVector GetWBoson(bool &status);
   std::vector<TLorentzVector*> GetNeutrinos(TLorentzVector* L, TLorentzVector* MET);
   virtual ~EventLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TString sampleName, TString ExpUncertaintyName);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

    
   Chi2_minimization* myMinimizer = new Chi2_minimization("MeV");
   NeutrinoBuilder* m_NeutrinoBuilder;
   
   TH1Fs *h_NBtags;
   TH1Fs *h_MET;
   TH1Fs *h_Lepton_Eta;
   TH1Fs *h_Lepton_Pt;
   TH1Fs *h_Mwt;

   TH1Fs *h_HT_bjets;
   TH1Fs *h_mVH;
   TH1Fs *h_DeltaPhi_HW;
   TH1Fs *h_DeltaEta_HW;
   TH1Fs *h_pTH;
   TH1Fs *h_pTH_over_mVH;
   TH1Fs *h_pTWplus;
   TH1Fs *h_pTW_over_mVH;
   TH1Fs *h_mH;
   TH1Fs *h_mWplus;
   TH1Fs *h_tagCategory;
   TH1Fs *h_mass_resolution;
   
   TH1Fs *h_NBtags_1b;
   TH1Fs *h_MET_1b;
   TH1Fs *h_Lepton_Eta_1b;
   TH1Fs *h_Lepton_Pt_1b;
   TH1Fs *h_Mwt_1b;

   TH1Fs *h_HT_bjets_1b;
   TH1Fs *h_mVH_1b;
   TH1Fs *h_DeltaPhi_HW_1b;
   TH1Fs *h_DeltaEta_HW_1b;
   TH1Fs *h_pTH_1b;
   TH1Fs *h_pTH_over_mVH_1b;
   TH1Fs *h_pTWplus_1b;
   TH1Fs *h_pTW_over_mVH_1b;
   TH1Fs *h_mH_1b;
   TH1Fs *h_mWplus_1b;
   TH1Fs *h_tagCategory_1b;
   TH1Fs *h_mass_resolution_1b; 

   TH1Fs *h_NBtags_1c;
   TH1Fs *h_MET_1c;
   TH1Fs *h_Lepton_Eta_1c;
   TH1Fs *h_Lepton_Pt_1c;
   TH1Fs *h_Mwt_1c;

   TH1Fs *h_HT_bjets_1c;
   TH1Fs *h_mVH_1c;
   TH1Fs *h_DeltaPhi_HW_1c;
   TH1Fs *h_DeltaEta_HW_1c;
   TH1Fs *h_pTH_1c;
   TH1Fs *h_pTH_over_mVH_1c;
   TH1Fs *h_pTWplus_1c;
   TH1Fs *h_pTW_over_mVH_1c;
   TH1Fs *h_mH_1c;
   TH1Fs *h_mWplus_1c;
   TH1Fs *h_tagCategory_1c;
   TH1Fs *h_mass_resolution_1c;

   TH1Fs *h_NBtags_1l;
   TH1Fs *h_MET_1l;
   TH1Fs *h_Lepton_Eta_1l;
   TH1Fs *h_Lepton_Pt_1l;
   TH1Fs *h_Mwt_1l;
   
   TH1Fs *h_HT_bjets_1l;
   TH1Fs *h_mVH_1l;
   TH1Fs *h_DeltaPhi_HW_1l;
   TH1Fs *h_DeltaEta_HW_1l;
   TH1Fs *h_pTH_1l;
   TH1Fs *h_pTH_over_mVH_1l;
   TH1Fs *h_pTWplus_1l;
   TH1Fs *h_pTW_over_mVH_1l;
   TH1Fs *h_mH_1l;
   TH1Fs *h_mWplus_1l;
   TH1Fs *h_tagCategory_1l;
   TH1Fs *h_mass_resolution_1l;

   vector<double>  m_EventWeights;
   vector<double>  m_ones;
   vector<TString> m_UncNames;
   vector<float> m_LepCh;
   vector<TString> mySel;

   std::vector<TLorentzVector> bQuarks;
   std::vector<TLorentzVector> LightQuarks;   
   std::vector<TLorentzVector> Leptons;
   std::vector<TLorentzVector> GenLevLeptons;
   std::vector<TLorentzVector> Jets;
   std::vector<TLorentzVector> FJets;
   std::vector<TLorentzVector> TrkJets;
   
   std::vector<int> Jets_PCbtag;
   std::vector<int> TrkJets_PCbtag;
   std::vector<int> nTaggedVRTrkJetsInFJet;

   TLorentzVector MET;
   TLorentzVector Higgs;
   TLorentzVector Wplus;
   TLorentzVector Higgs_tr;
   TLorentzVector Wminus;
   TLorentzVector Higgs_LV;
   TLorentzVector Wplus_LV;
   TLorentzVector Higgs_p1;
   TLorentzVector Higgs_p2;
   TLorentzVector Wplus_p1;
   TLorentzVector Wplus_p2;
   bool   found_Wplus_had;
   bool   found_Wplus_lep;
   bool   found_Higgs;
   int    m_NTags;

   int    m_NTags_caloJ;
   int    m_NTags_trkJ;
   int    m_NTags_Higgs;
   int    m_ntagsOutside;
   int    m_bTagCategory;
   int    m_btagCategoryBin;

   float m_mWT;
   float m_Lep_pT;
   float m_Lep_Eta;
   float m_mVH;
   float m_DeltaPhi_HW;
   float m_DeltaEta_HW;
   float m_MassTruth;
   float m_MET;
   float m_H_mass;    
   float m_H_pT;
   float m_Wp_mass;
   float m_Wp_pT;
   int m_flag_tt;
   float m_weight_scale;

   
};

#endif

#ifdef EventLoop_cxx
EventLoop::EventLoop(TTree *tree, TString sampleName, TString MCDataPeriode, TString ExpUncertaintyName, TString WP) : fChain(0) {

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
	std::cerr << "Error in EventLoop::EventLoop(): tree is nullptr" << std::endl;
	return;
   }

   if(WP == "85p"){
      m_btagCategoryBin        = 2;
   }
   if(WP == "77p"){
      m_btagCategoryBin        = 3;
   } 
   if(WP == "70p"){
      m_btagCategoryBin        = 4;
   }
   if(WP == "60p"){
      m_btagCategoryBin        = 5; 
   }

   if(sampleName.Contains("tt_PP8.root") || sampleName.Contains("tt_PP8filtered.root"))
   {
     m_flag_tt = 1;
   }
   else
   {
    m_flag_tt = 0;
   }

   if(sampleName=="hp400_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/40524.20;
   }

   if(sampleName=="hp800_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/2126.78;
   }

   if(sampleName=="hp1200_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/304.93;
   }

   if(sampleName=="hp1600_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/46.37;
   }

   if(sampleName=="hp2000_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/14.89;
   }

   if(sampleName=="hp2500_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/3.20;
   }

   if(sampleName=="hp3000_AFII.root" && MCDataPeriode == "MC16a"){
     m_weight_scale = 1/0.512;
   }

   

   if(sampleName=="hp400_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/30018.2;
   }

   if(sampleName=="hp800_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/3178.54;
   }

   if(sampleName=="hp1200_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/199.62;
   }

   if(sampleName=="hp1600_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/71.57;
   }

   if(sampleName=="hp2000_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/19.17;
   }

   if(sampleName=="hp2500_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/2.62;
   }

   if(sampleName=="hp3000_AFII.root" && MCDataPeriode == "MC16d"){
     m_weight_scale = 1/0.46;
   }


   if(sampleName=="hp400_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/51948;
   }

   if(sampleName=="hp800_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/3633.95;
   }

   if(sampleName=="hp1200_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/311.36;
   }

   if(sampleName=="hp1600_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/80.22;
   }

   if(sampleName=="hp2000_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/27.57;
   }

   if(sampleName=="hp2500_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/3.35;
   }

   if(sampleName=="hp3000_AFII.root" && MCDataPeriode == "MC16e"){
     m_weight_scale = 1/1.25;
   }


   Init(tree, sampleName, ExpUncertaintyName);
}

EventLoop::~EventLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventLoop::Init(TTree *tree, TString sampleName, TString ExpUncertaintyName){
   

   m_NeutrinoBuilder = new NeutrinoBuilder("MeV");
   sampleName.Resize(sampleName.Index(".root"));

   m_UncNames   = {""};   
   
   mySel        = {"Merged_SR"};

   h_MET                    = new TH1Fs(sampleName+"_MET", "",        30, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Eta             = new TH1Fs(sampleName+"_Lepton_Eta", "", 25, -2.5, 2.5, mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Pt              = new TH1Fs(sampleName+"_Lepton_Pt", "",  25, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_NBtags                 = new TH1Fs(sampleName+"_nBTags", "",     7, -0.5, 6.5,  mySel, m_UncNames, ExpUncertaintyName);
   h_Mwt                    = new TH1Fs(sampleName+"_Mwt", "",        30,  0, 300,  mySel, m_UncNames, ExpUncertaintyName);
   h_mVH                    = new TH1Fs(sampleName+"_mVH", "",       35, 0, 3500,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH                    = new TH1Fs(sampleName+"_pTH", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_over_mVH           = new TH1Fs(sampleName+"_pTH_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTWplus                = new TH1Fs(sampleName+"_pTWplus", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTW_over_mVH      	    = new TH1Fs(sampleName+"_pTW_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_mH                     = new TH1Fs(sampleName+"_mH", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName); 
   h_mWplus                 = new TH1Fs(sampleName+"_mWplus", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName);
   
   
   h_MET_1b                   = new TH1Fs(sampleName+"_1B_MET", "",        30, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Eta_1b             = new TH1Fs(sampleName+"_1B_Lepton_Eta", "", 25, -2.5, 2.5, mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Pt_1b             = new TH1Fs(sampleName+"_1B_Lepton_Pt", "",  25, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_NBtags_1b                = new TH1Fs(sampleName+"_1B_nBTags", "",     7, -0.5, 6.5,  mySel, m_UncNames, ExpUncertaintyName);
   h_Mwt_1b                    = new TH1Fs(sampleName+"_1B_Mwt", "",        30,  0, 300,  mySel, m_UncNames, ExpUncertaintyName);
   h_mVH_1b                    = new TH1Fs(sampleName+"_1B_mVH", "",       35, 0, 3500,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_1b                    = new TH1Fs(sampleName+"_1B_pTH", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_over_mVH_1b          = new TH1Fs(sampleName+"_1B_pTH_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTWplus_1b               = new TH1Fs(sampleName+"_1B_pTWplus", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTW_over_mVH_1b      	    = new TH1Fs(sampleName+"_1B_pTW_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_mH_1b                     = new TH1Fs(sampleName+"_1B_mH", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName); 
   h_mWplus_1b                 = new TH1Fs(sampleName+"_1B_mWplus", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName);
   
   
   h_MET_1c                   = new TH1Fs(sampleName+"_1C_MET", "",        30, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Eta_1c             = new TH1Fs(sampleName+"_1C_Lepton_Eta", "", 25, -2.5, 2.5, mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Pt_1c             = new TH1Fs(sampleName+"_1C_Lepton_Pt", "",  25, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_NBtags_1c                 = new TH1Fs(sampleName+"_1C_nBTags", "",     7, -0.5, 6.5,  mySel, m_UncNames, ExpUncertaintyName);
   h_Mwt_1c                    = new TH1Fs(sampleName+"_1C_Mwt", "",        30,  0, 300,  mySel, m_UncNames, ExpUncertaintyName);
   h_mVH_1c                    = new TH1Fs(sampleName+"_1C_mVH", "",       35, 0, 3500,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_1c                    = new TH1Fs(sampleName+"_1C_pTH", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_over_mVH_1c           = new TH1Fs(sampleName+"_1C_pTH_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTWplus_1c                = new TH1Fs(sampleName+"_1C_pTWplus", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTW_over_mVH_1c     	    = new TH1Fs(sampleName+"_1C_pTW_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_mH_1c                     = new TH1Fs(sampleName+"_1C_mH", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName); 
   h_mWplus_1c                 = new TH1Fs(sampleName+"_1C_mWplus", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName);
   

   h_MET_1l                    = new TH1Fs(sampleName+"_1L_MET", "",        30, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Eta_1l             = new TH1Fs(sampleName+"_1L_Lepton_Eta", "", 25, -2.5, 2.5, mySel, m_UncNames, ExpUncertaintyName);
   h_Lepton_Pt_1l             = new TH1Fs(sampleName+"_1L_Lepton_Pt", "",  25, 0, 600,    mySel, m_UncNames, ExpUncertaintyName);
   h_NBtags_1l                 = new TH1Fs(sampleName+"_1L_nBTags", "",     7, -0.5, 6.5,  mySel, m_UncNames, ExpUncertaintyName);
   h_Mwt_1l                    = new TH1Fs(sampleName+"_1L_Mwt", "",        30,  0, 300,  mySel, m_UncNames, ExpUncertaintyName);
   h_mVH_1l                   = new TH1Fs(sampleName+"_1L_mVH", "",       35, 0, 3500,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_1l                   = new TH1Fs(sampleName+"_1L_pTH", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTH_over_mVH_1l          = new TH1Fs(sampleName+"_1L_pTH_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTWplus_1l                = new TH1Fs(sampleName+"_1L_pTWplus", "",25, 0, 1000,    mySel, m_UncNames, ExpUncertaintyName);
   h_pTW_over_mVH_1l     	    = new TH1Fs(sampleName+"_1L_pTW_over_mVH", "",25, 0, 1,    mySel, m_UncNames, ExpUncertaintyName);
   h_mH_1l                     = new TH1Fs(sampleName+"_1L_mH", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName); 
   h_mWplus_1l                 = new TH1Fs(sampleName+"_1L_mWplus", "",30, 50, 200,    mySel, m_UncNames, ExpUncertaintyName);
   
   

   // Set object pointer
   el_pt                 = 0;
   el_eta                = 0;
   el_phi                = 0;
   el_e                  = 0;
   el_charge             = 0;

   mu_pt                 = 0;
   mu_eta                = 0;
   mu_phi                = 0;
   mu_e                  = 0;
   mu_charge             = 0;

   signal_Jet_E          = 0;
   signal_Jet_PT         = 0;
   signal_Jet_Eta        = 0;
   signal_Jet_Phi        = 0;
   signal_Jet_truthflav  = 0; 
   signal_Jet_tagWeightBin_DL1r_Continuous = 0;

   tjet_e          = 0;
   tjet_pt         = 0;
   tjet_eta        = 0;
   tjet_phi        = 0;
   tjet_tagWeightBin_DL1r_Continuous = 0;

   ljet_m          = 0;
   ljet_pt         = 0;
   ljet_eta        = 0;
   ljet_phi        = 0;
   ljet_truthLabel = 0;

   ljet_D2                   = 0;
   ljet_C2                   = 0;
   ljet_Xbb2020v3_Higgs      = 0;
   ljet_Xbb2020v3_QCD        = 0;
   ljet_Xbb2020v3_Top        = 0;
   
   truth_pt              = 0;
   truth_eta             = 0;
   truth_phi             = 0;
   truth_m               = 0;
   truth_pdgid           = 0; 
   truth_status          = 0;
   truth_tthbb_info      = 0; 

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("nJets", &nJet, &b_nJet);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("nPrimaryVtx", &nPrimaryVtx, &b_nPrimaryVtx);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("randomRunNumber", &randomRunNumber, &b_randomRunNumber);
   
   
   fChain->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
   fChain->SetBranchAddress("weight_jvt", &weight_jvt, &b_weight_jvt);
   fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   fChain->SetBranchAddress("weight_leptonSF", &weight_leptonSF, &b_weight_leptonSF);
   fChain->SetBranchAddress("weight_trackjet_bTagSF_DL1r_Continuous", &weight_trackjet_bTagSF_DL1r_Continuous, &b_weight_trackjet_bTagSF_DL1r_Continuous);
    
   fChain->SetBranchAddress("nTaus",      &nTaus,      &b_nTaus);
   fChain->SetBranchAddress("met_met", &met, &b_MET);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_MET_phi);

   
   fChain->SetBranchAddress("el_pt",   &el_pt,  &b_el_pt);
   fChain->SetBranchAddress("el_eta",  &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi",  &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_e",    &el_e,   &b_el_e);
   fChain->SetBranchAddress("el_charge",        &el_charge,       &b_el_charge);

   fChain->SetBranchAddress("mu_pt",   &mu_pt,  &b_mu_pt);
   fChain->SetBranchAddress("mu_eta",  &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi",  &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_e",    &mu_e,   &b_mu_e);
   fChain->SetBranchAddress("mu_charge",    &mu_charge,   &b_mu_charge);

   fChain->SetBranchAddress("ljet_m", &ljet_m, &b_ljet_m);
   fChain->SetBranchAddress("ljet_pt", &ljet_pt, &b_ljet_pt);
   fChain->SetBranchAddress("ljet_eta", &ljet_eta, &b_ljet_eta);
   fChain->SetBranchAddress("ljet_phi", &ljet_phi, &b_ljet_phi);
   fChain->SetBranchAddress("ljet_truthLabel", &ljet_truthLabel, &b_ljet_truthLabel);

   fChain->SetBranchAddress("ljet_D2", &ljet_D2, &b_ljet_D2);
   fChain->SetBranchAddress("ljet_C2", &ljet_C2, &b_ljet_C2);
   fChain->SetBranchAddress("ljet_Xbb2020v3_Higgs", &ljet_Xbb2020v3_Higgs, &b_ljet_Xbb2020v3_Higgs);
   fChain->SetBranchAddress("ljet_Xbb2020v3_QCD", &ljet_Xbb2020v3_QCD, &b_ljet_Xbb2020v3_QCD);
   fChain->SetBranchAddress("ljet_Xbb2020v3_Top", &ljet_Xbb2020v3_Top, &b_ljet_Xbb2020v3_Top);
   
   fChain->SetBranchAddress("tjet_e",   &tjet_e, &b_tjet_e);
   fChain->SetBranchAddress("tjet_pt",  &tjet_pt, &b_tjet_pt);
   fChain->SetBranchAddress("tjet_eta", &tjet_eta, &b_tjet_eta);
   fChain->SetBranchAddress("tjet_phi", &tjet_phi, &b_tjet_phi);
   fChain->SetBranchAddress("tjet_tagWeightBin_DL1r_Continuous", &tjet_tagWeightBin_DL1r_Continuous, &b_tjet_tagWeightBin_DL1r_Continuous);
   
   fChain->SetBranchAddress("jet_e",   &signal_Jet_E, &b_signal_Jet_E);
   fChain->SetBranchAddress("jet_pt",  &signal_Jet_PT, &b_signal_Jet_PT);
   fChain->SetBranchAddress("jet_eta", &signal_Jet_Eta, &b_signal_Jet_Eta);
   fChain->SetBranchAddress("jet_phi", &signal_Jet_Phi, &b_signal_Jet_Phi);
   fChain->SetBranchAddress("jet_truthflav", &signal_Jet_truthflav, &b_signal_Jet_truthflav); //only for MC
   fChain->SetBranchAddress("jet_tagWeightBin_DL1r_Continuous", &signal_Jet_tagWeightBin_DL1r_Continuous, &b_signal_Jet_tagWeightBin_DL1r_Continuous);
   

   fChain->SetBranchAddress("truth_pt", &truth_pt, &b_truth_pt);
   fChain->SetBranchAddress("truth_eta", &truth_eta, &b_truth_eta);
   fChain->SetBranchAddress("truth_phi", &truth_phi, &b_truth_phi);
   fChain->SetBranchAddress("truth_m", &truth_m, &b_truth_m);
   fChain->SetBranchAddress("truth_pdgid", &truth_pdgid, &b_truth_pdgid);
   fChain->SetBranchAddress("truth_status", &truth_status, &b_truth_status);
   fChain->SetBranchAddress("truth_tthbb_info", &truth_tthbb_info, &b_truth_tthbb_info); 
   fChain->SetBranchAddress("HF_SimpleClassification", &HF_SimpleClassification, &b_HF_SimpleClassification);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("TopHeavyFlavorFilterFlag", &TopHeavyFlavorFilterFlag, &b_TopHeavyFlavorFilterFlag);
   
   Notify();
}

Bool_t EventLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventLoop_cxx

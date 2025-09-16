#define EventLoop_cxx
#define _USE_MATH_DEFINES
#include <cmath>
#include "EventLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include <TROOT.h>
#include <TRint.h>

void EventLoop::Loop(){
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000000 == 0)std::cout<<"Processing "<<jentry<<" events!!!"<<std::endl;
       
      for(auto sel : mySel){
          pass_sel[sel] = false;
      }

      m_NTags         = 0;
      m_NTags_caloJ    = 0;
      m_NTags_trkJ     = 0;
      m_NTags_Higgs    = 0;
      m_ntagsOutside   = 0;
      m_MassTruth      = 0;

      EventWeight      = m_weight_scale*weight_pileup*weight_leptonSF*weight_jvt*weight_mc*weight_trackjet_bTagSF_DL1r_Continuous*(randomRunNumber<=311481 ? 36207.66 : (randomRunNumber<=340453 ? 44307.4 : 58450.1));
      
      m_EventWeights.clear();
      m_EventWeights.push_back(EventWeight);
      
      m_ones.clear();
      m_ones.push_back(1.0);
      MET.SetPtEtaPhiM(met,0.,met_phi,0.);
      SetJetVectors(); 
      SetLeptonVectors();
      SetTruthParticles();

      bool status_W = false;
      m_mWT         = GetMwt();
      m_MassTruth   = GetTruthMass();
      Wminus        = GetWBoson(status_W); //Wminus from top in association with charged Higgs, built from lepton and neutrino
      bool passed_merged_preselection = PassEventSelectionBoosted();
      if(!passed_merged_preselection)continue;
      
        if(passed_merged_preselection == true)
    {
                   pass_sel["Merged_SR"]  = true;
	 } 
   
      if(pass_sel["Merged_SR"]==true)
	   {  
         if(m_flag_tt == 1 && (((mcChannelNumber==410470) && (TopHeavyFlavorFilterFlag==0)) || ((mcChannelNumber==411073 || mcChannelNumber==411076) && (TopHeavyFlavorFilterFlag==1)) || ((mcChannelNumber==411074 || mcChannelNumber==411077) && (TopHeavyFlavorFilterFlag==2)) || ((mcChannelNumber==411075 || mcChannelNumber==411078) && (TopHeavyFlavorFilterFlag==3)))) {      
            if(HF_SimpleClassification == 1){     
         h_MET_1b->Fill(MET.Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Eta_1b->Fill(Leptons.at(0).Eta(),m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Pt_1b->Fill(Leptons.at(0).Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_NBtags_1b->Fill(m_NTags,m_EventWeights, pass_sel, m_NTags);
         h_Mwt_1b->Fill(m_mWT,m_EventWeights, pass_sel, m_NTags);
         h_mVH_1b->Fill(m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTH_1b->Fill(Higgs.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTWplus_1b->Fill(Wplus.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTH_over_mVH_1b->Fill(Higgs.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTW_over_mVH_1b->Fill(Wplus.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_mH_1b->Fill(Higgs.M()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_mWplus_1b->Fill(Wplus.M()*0.001, m_EventWeights, pass_sel, m_NTags);
            }

            else if(HF_SimpleClassification == -1){
         h_MET_1c->Fill(MET.Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Eta_1c->Fill(Leptons.at(0).Eta(),m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Pt_1c->Fill(Leptons.at(0).Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_NBtags_1c->Fill(m_NTags,m_EventWeights, pass_sel, m_NTags);
         h_Mwt_1c->Fill(m_mWT,m_EventWeights, pass_sel, m_NTags);
         h_mVH_1c->Fill(m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTH_1c->Fill(Higgs.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTWplus_1c->Fill(Wplus.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTH_over_mVH_1c->Fill(Higgs.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTW_over_mVH_1c->Fill(Wplus.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_mH_1c->Fill(Higgs.M()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_mWplus_1c->Fill(Wplus.M()*0.001, m_EventWeights, pass_sel, m_NTags);
            }

            else if(HF_SimpleClassification == 0){
         h_MET_1l->Fill(MET.Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Eta_1l->Fill(Leptons.at(0).Eta(),m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Pt_1l->Fill(Leptons.at(0).Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_NBtags_1l->Fill(m_NTags,m_EventWeights, pass_sel, m_NTags);
         h_Mwt_1l->Fill(m_mWT,m_EventWeights, pass_sel, m_NTags);
         h_mVH_1l->Fill(m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTH_1l->Fill(Higgs.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTWplus_1l->Fill(Wplus.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTH_over_mVH_1l->Fill(Higgs.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTW_over_mVH_1l->Fill(Wplus.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_mH_1l->Fill(Higgs.M()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_mWplus_1l->Fill(Wplus.M()*0.001, m_EventWeights, pass_sel, m_NTags);
            }
         } 

	 else{

         h_MET->Fill(MET.Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Eta->Fill(Leptons.at(0).Eta(),m_EventWeights, pass_sel, m_NTags);
         h_Lepton_Pt->Fill(Leptons.at(0).Pt()*0.001,m_EventWeights, pass_sel, m_NTags);
         h_NBtags->Fill(m_NTags,m_EventWeights, pass_sel, m_NTags);
         h_Mwt->Fill(m_mWT,m_EventWeights, pass_sel, m_NTags);
         h_mVH->Fill(m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTH->Fill(Higgs.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTWplus->Fill(Wplus.Pt()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_pTH_over_mVH->Fill(Higgs.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_pTW_over_mVH->Fill(Wplus.Pt()*0.001/m_mVH, m_EventWeights, pass_sel, m_NTags);
         h_mH->Fill(Higgs.M()*0.001, m_EventWeights, pass_sel, m_NTags);
         h_mWplus->Fill(Wplus.M()*0.001, m_EventWeights, pass_sel, m_NTags);
            }
      } 
   } 
}


//********Use the flags directly from tthbb info**************

void EventLoop::SetTruthParticles(){
   TLorentzVector truthParticle;
   bool contains_truth_Higgs = false;
   bool contains_truth_W     = false;
   int  n_attempts           = 0;
   bool bfromH               = false;
   bool bbarfromH            = false;
   bool finst                = false;
   bool fromt                = false;
   bool fromtbar             = false;
   bool tdirec               = false;
   bool assochp              = false;
   bool fhplus               = false;
   found_Wplus_had           = false;
   found_Wplus_lep           = false;
   found_Higgs               = false;

   while(!found_Higgs || !(found_Wplus_lep || found_Wplus_had)){
      bQuarks.clear();
      LightQuarks.clear();
      GenLevLeptons.clear();

      for(unsigned int i = 0; i < truth_pt->size(); i++){
         
         truthParticle.SetPtEtaPhiM(truth_pt->at(i),truth_eta->at(i),truth_phi->at(i),truth_m->at(i));
         //flags based on final state decay of particles
         bfromH = truth_tthbb_info->at(i)&(1<<0);  //its a b from Higgs
         bbarfromH = truth_tthbb_info->at(i)&(1<<1); //its a bbar from Higgs
         finst = truth_tthbb_info->at(i)&( (Long64_t) 1<<50); //final state flag
         fromt = truth_tthbb_info->at(i)&(1<<28); // decay product from top quark
         fromtbar = truth_tthbb_info->at(i)&(1<<29); // decay product from anti-top quark
         assochp = truth_tthbb_info->at(i)&(1<<8); // particle associated with charged Higgs
         tdirec = truth_tthbb_info->at(i)&(1<<17); //direct decay from top quark
         
         if(fabs(truth_pdgid->at(i)) <= 4 && fabs(truth_pdgid->at(i)) != 5 && finst == true && fromt == false && fromtbar == false && assochp == false){   
             LightQuarks.push_back(truthParticle);
         }
         
         if( fabs(truth_pdgid->at(i)) <= 16 && fabs(truth_pdgid->at(i)) >= 11 && fabs(truth_pdgid->at(i)) != 5 && finst == true && fromt == false && fromtbar == false && assochp == false){
             GenLevLeptons.push_back(truthParticle);
         }
            
         if(truth_pdgid->at(i) == 24){
             Wplus.SetPtEtaPhiM(truth_pt->at(i),truth_eta->at(i),truth_phi->at(i),truth_m->at(i));
             contains_truth_W     = true;
         }
      }

      for(unsigned int i = 0; i < truth_pt->size(); i++){

          if(truth_pdgid->at(i) == 25){
             Higgs_tr.SetPtEtaPhiM(truth_pt->at(i),truth_eta->at(i),truth_phi->at(i),truth_m->at(i));
             found_Higgs  = true;
         }

      }

      
      for(unsigned int i = 0; i < LightQuarks.size(); i++){
         float dRmin = 42.0; 
         for(unsigned int j = 0; j < LightQuarks.size(); j++){   
           if(i == j)continue;
            //if(contains_truth_W ){       
                Wplus_p1         = LightQuarks.at(i);
                Wplus_p2         = LightQuarks.at(j);
                found_Wplus_had  = true;
                //std::cout<<"Found Hadronic Wplus"<<std::endl;
           //}
         }
      }
      
     
      for(unsigned int i = 0; i < GenLevLeptons.size(); i++){
         float dRmin = 42.0;
          for(unsigned int j = 0; j < GenLevLeptons.size(); j++){  
           if(i == j)continue;
           //if(contains_truth_W){
                Wplus_p1         = GenLevLeptons.at(i);
                Wplus_p2         = GenLevLeptons.at(j);
                found_Wplus_lep  = true;
                //std::cout<<"Found Leptonic Wplus"<<std::endl;
           //}
         }
      } 

      if(!found_Higgs){
         if(n_attempts == 2)break;
         n_attempts++;
      }
      
      if(!(found_Wplus_lep || found_Wplus_had)){
         if(n_attempts == 2)break;
         n_attempts++;
      } 
   }
}



int EventLoop::GetBTagCategoryShort(int NTags_InHiggsJet, int NTags_OutsideHiggsJet){
     int category = -1;
     if(NTags_InHiggsJet == 0 && NTags_OutsideHiggsJet == 2)category=2;
     if(NTags_InHiggsJet == 1 && NTags_OutsideHiggsJet == 1)category=2;
     if(NTags_InHiggsJet == 2 && NTags_OutsideHiggsJet == 0)category=2;
     if(NTags_InHiggsJet == 0 && NTags_OutsideHiggsJet == 3)category=3;
     if(NTags_InHiggsJet == 1 && NTags_OutsideHiggsJet == 2)category=3;
     if(NTags_InHiggsJet == 2 && NTags_OutsideHiggsJet == 1)category=3;
     if(NTags_InHiggsJet == 0 && NTags_OutsideHiggsJet >= 4)category=4;
     if(NTags_InHiggsJet == 1 && NTags_OutsideHiggsJet >= 3)category=4;
     if(NTags_InHiggsJet >= 2 && NTags_OutsideHiggsJet >= 2)category=4;
     return category;
}

double EventLoop::GetTruthMass(){
     if(found_Higgs && (found_Wplus_had || found_Wplus_lep)){
         return (Higgs_tr+(Wplus_p1+Wplus_p2)).M()*0.001; //check the naming: maybe Higgs used at reco level also, if thats case use Higgs_tr at truth level
     }
     return 999;
}

//Function to find the Fat Jet pair for Higgs and W
bool EventLoop::FindFJetPair(){
      bool status   = false;
      Higgs  = FJets.at(0).M() > FJets.at(1).M() ? FJets.at(0) : FJets.at(1);
      Wplus  = FJets.at(0).M() < FJets.at(1).M() ? FJets.at(0) : FJets.at(1);
      m_NTags_Higgs = FJets.at(0).M() > FJets.at(1).M() ? nTaggedVRTrkJetsInFJet.at(0) : nTaggedVRTrkJetsInFJet.at(1);
      m_ntagsOutside     = m_NTags_trkJ - (nTaggedVRTrkJetsInFJet.at(0) + nTaggedVRTrkJetsInFJet.at(1));
      m_NTags            = GetBTagCategoryShort(m_NTags_Higgs,m_ntagsOutside);
      if(Higgs.Pt() > 250000)status = true;
      m_DeltaPhi_HW      = fabs(Wplus.DeltaPhi(Higgs));
      m_mVH              = (Wplus+Higgs).M()*0.001;
      return status;
}

TLorentzVector EventLoop::GetWBoson(bool &status){
   status = false;
   TLorentzVector Wleptonic;
   if(Leptons.size() == 0) return Wleptonic;
   std::vector<TLorentzVector*> neutrinoVector = GetNeutrinos(&Leptons.at(0),&MET);
   for(auto neutrino : neutrinoVector){
       Wleptonic = (*neutrino + Leptons.at(0));
       status = true;
   }
   return Wleptonic;
}

std::vector<TLorentzVector*> EventLoop::GetNeutrinos(TLorentzVector* L, TLorentzVector* MET){
   std::vector<TLorentzVector*> neutrinoVector;
   neutrinoVector      = m_NeutrinoBuilder->candidatesFromWMass_Rotation(L, MET, true);
   bool m_isRotatedSol = m_NeutrinoBuilder->m_isRotated;
   double m_r          = m_NeutrinoBuilder->m_r;
   return neutrinoVector;
}

//Function to select the boosted events
bool EventLoop::PassEventSelectionBoosted(){
    if(MET.Pt() < 30000)return false;  
    if(Leptons.size() != 1)return false;
    if(FJets.size() < 2)return false;
    if(FJets.at(0).Pt() < 200000)return false;
    if(FJets.at(1).Pt() < 200000)return false;
    if(FJets.at(0).DeltaR(Leptons.at(0)) < 1.0)return false;
    if(FJets.at(1).DeltaR(Leptons.at(0)) < 1.0)return false;
    if(!FindFJetPair())return false;
    if(m_DeltaPhi_HW < 2.5)return false;
    return true;
}


double EventLoop::GetMwt(){
	return Leptons.size() > 0 ? (0.001*sqrt(2. * Leptons.at(0).Pt() * MET.Pt() * (1. - cos(Leptons.at(0).DeltaPhi(MET))))) : 0.;
}


void EventLoop::SetLeptonVectors(){
   Leptons.clear();
   m_LepCh.clear();
   TLorentzVector lepton;
   bool trigger_matched = false;
   for(int i=0; i < mu_pt->size(); i++) {
        if(mu_pt->at(i)*0.001 < 27.)continue;
        lepton.SetPtEtaPhiE(mu_pt->at(i), mu_eta->at(i), mu_phi->at(i), mu_e->at(i));
        Leptons.push_back(lepton);
        m_LepCh.push_back(mu_charge->at(i));
   }
   for(int i=0; i < el_pt->size(); i++) {
        if(el_pt->at(i)*0.001 < 27.)continue;
        lepton.SetPtEtaPhiE(el_pt->at(i), el_eta->at(i), el_phi->at(i), el_e->at(i));
        Leptons.push_back(lepton);
        m_LepCh.push_back(el_charge->at(i));
   }
}

void EventLoop::Calc_dR_jl(float &dR_jl_min, float &dR_jl_max){
     dR_jl_min = 42.0;
     dR_jl_max =  0.0;
     for(int i=0; i < (Jets.size() <= 6 ? Jets.size() : 6); i++) {
          if(dR_jl_min > Jets.at(i).DeltaR(Leptons.at(0))){
              dR_jl_min = Jets.at(i).DeltaR(Leptons.at(0));
          }
       	  if(dR_jl_max < Jets.at(i).DeltaR(Leptons.at(0))){
       	      dR_jl_max	= Jets.at(i).DeltaR(Leptons.at(0));
       	  }
     }
}

void EventLoop::SetJetVectors(){
    TrkJets.clear();
    TrkJets_PCbtag.clear();
    FJets.clear();
    nTaggedVRTrkJetsInFJet.clear();
    TLorentzVector jet, trkjet;
    for(int i = 0; i < ljet_m->size(); i++) {
        if(ljet_pt->at(i) < 200000)continue;
        if(ljet_m->at(i)  <  50000)continue;
	jet.SetPtEtaPhiM(ljet_pt->at(i), ljet_eta->at(i), ljet_phi->at(i), ljet_m->at(i));
	FJets.push_back(jet);
        int nMatchedAndTaggedJets = 0;        
        for(int i=0; i < tjet_pt->size(); i++) {
             trkjet.SetPtEtaPhiE(tjet_pt->at(i), tjet_eta->at(i), tjet_phi->at(i), tjet_e->at(i));             
             if(jet.DeltaR(trkjet) < 1.0 && tjet_tagWeightBin_DL1r_Continuous->at(i) >= m_btagCategoryBin) nMatchedAndTaggedJets++;
        }
        nTaggedVRTrkJetsInFJet.push_back(nMatchedAndTaggedJets); 
    }
    for(int i=0; i < tjet_pt->size(); i++) {
        jet.SetPtEtaPhiE(tjet_pt->at(i), tjet_eta->at(i), tjet_phi->at(i), tjet_e->at(i));
        TrkJets.push_back(jet);
        if(tjet_tagWeightBin_DL1r_Continuous->at(i) >= m_btagCategoryBin) {
            TrkJets_PCbtag.push_back(1);
            m_NTags_trkJ++;
	}else{
            TrkJets_PCbtag.push_back(0);
        }
    }
    Sort_Jets(&FJets,&nTaggedVRTrkJetsInFJet);
}


void EventLoop :: Sort_Jets(std::vector<TLorentzVector> *Jets, std::vector<int> *is_tagged){
    bool bDone = false;
    while (!bDone){
     	bDone = true;
        for (unsigned int i = 0; i < Jets->size()-1 && Jets->size()>0; i++){
            if ( Jets->at(i).Pt() < Jets->at(i+1).Pt()){
             	TLorentzVector tmp_jet   = Jets->at(i);
                Jets->at(i)              = Jets->at(i+1);
                Jets->at(i+1)            = tmp_jet;
                int tmp_tagging_info     = is_tagged->at(i);
                is_tagged->at(i)         = is_tagged->at(i+1);
                is_tagged->at(i+1)       = tmp_tagging_info;
                bDone = false;
            }
	}
    }
}

void EventLoop::Write(TDirectory *dir, std::string dirname) {

    h_MET->Write(dir, ("MET"));
    h_Mwt->Write(dir, ("Mwt"));
    h_Lepton_Eta->Write(dir, ("Lepton_Eta"));
    h_Lepton_Pt->Write(dir, ("Lepton_Pt"));
    h_NBtags->Write(dir, ("nBTags"));
    h_mVH->Write(dir, ("mVH"));
    h_pTH->Write(dir, ("pTH"));
    h_pTH_over_mVH->Write(dir, ("pTH_over_mVH"));
    h_pTWplus->Write(dir, ("pTWplus"));
    h_pTW_over_mVH->Write(dir, ("pTW_over_mVH"));
    h_mH->Write(dir, ("mH"));
    h_mWplus->Write(dir, ("mWplus"));
    h_tagCategory->Write(dir, ("BtagCategory"));
    
    h_MET_1b->Write(dir, ("MET_1b"));
    h_Mwt_1b->Write(dir, ("Mwt_1b"));
    h_Lepton_Eta_1b->Write(dir, ("Lepton_Eta_1b"));
    h_Lepton_Pt_1b->Write(dir, ("Lepton_Pt_1b"));
    h_NBtags_1b->Write(dir, ("nBTags_1b"));
    h_mVH_1b->Write(dir, ("mVH_1b"));
    h_pTH_1b->Write(dir, ("pTH_1b"));
    h_pTH_over_mVH_1b->Write(dir, ("pTH_over_mVH_1b"));
    h_pTWplus_1b->Write(dir, ("pTWplus_1b"));
    h_pTW_over_mVH_1b->Write(dir, ("pTW_over_mVH_1b"));
    h_mH_1b->Write(dir, ("mH_1b"));
    h_mWplus_1b->Write(dir, ("mWplus_1b"));
    h_tagCategory_1b->Write(dir, ("BtagCategory_1b"));

    h_MET_1c->Write(dir, ("MET_1c"));
    h_Mwt_1c->Write(dir, ("Mwt_1c"));
    h_Lepton_Eta_1c->Write(dir, ("Lepton_Eta_1c"));
    h_Lepton_Pt_1c->Write(dir, ("Lepton_Pt_1c"));
    h_NBtags_1c->Write(dir, ("nBTags_1c"));
    h_mVH_1c->Write(dir, ("mVH_1c"));
    h_pTH_1c->Write(dir, ("pTH_1c"));
    h_pTH_over_mVH_1c->Write(dir, ("pTH_over_mVH_1c"));
    h_pTWplus_1c->Write(dir, ("pTWplus_1c"));
    h_pTW_over_mVH_1c->Write(dir, ("pTW_over_mVH_1c"));
    h_mH_1c->Write(dir, ("mH_1c"));
    h_mWplus_1c->Write(dir, ("mWplus_1c"));
    h_tagCategory_1c->Write(dir, ("BtagCategory_1c"));
    
    h_MET_1l->Write(dir, ("MET_1l"));
    h_Mwt_1l->Write(dir, ("Mwt_1l"));
    h_Lepton_Eta_1l->Write(dir, ("Lepton_Eta_1l"));
    h_Lepton_Pt_1l->Write(dir, ("Lepton_Pt_1l"));
    h_NBtags_1l->Write(dir, ("nBTags_1l"));
    h_mVH_1l->Write(dir, ("mVH_1l"));
    h_pTH_1l->Write(dir, ("pTH_1l"));
    h_pTH_over_mVH_1l->Write(dir, ("pTH_over_mVH_1l"));
    h_pTWplus_1l->Write(dir, ("pTWplus_1l"));
    h_pTW_over_mVH_1l->Write(dir, ("pTW_over_mVH_1l"));
    h_mH_1l->Write(dir, ("mH_1l"));
    h_mWplus_1l->Write(dir, ("mWplus_1l"));
    h_tagCategory_1l->Write(dir, ("BtagCategory_1l"));
}

void EventLoop::WriteTreesToFile(TFile *outFile) {
    outFile->cd();
}

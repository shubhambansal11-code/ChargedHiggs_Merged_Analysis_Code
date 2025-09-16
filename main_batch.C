#include "main/EventLoop.h"
#include "TApplication.h"
#include <TSystemDirectory.h>
#include <TList.h>
#include <TCollection.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

int main(int argc, char** argv){
  clock_t begin = clock();

  TString path          =  TString(argv[1]);
  TString SampleName    =  TString(argv[2]);
  TString WP            =  TString(argv[3]);
  TString OUTPUTDIR     =  TString(argv[4]);
  
  TString MCDataPeriode = "";
  std::vector<std::string> TreeNames;
  TreeNames.push_back("nominal_Loose");
  TreeNames.push_back("sumWeights");
  if(path.Contains("MC16a"))MCDataPeriode = "MC16a";
  if(path.Contains("MC16d"))MCDataPeriode = "MC16d";
  if(path.Contains("MC16e"))MCDataPeriode = "MC16e";
  
  TString OutFile = SampleName+"_"+WP+"_"+MCDataPeriode+".root";
  
  TString OutFileName = OUTPUTDIR + "/"+ OutFile;
  
  TFile* outfile = TFile::Open(OutFileName,"RECREATE");
  for(unsigned int i=0; i < TreeNames.size(); i++){
     
     std::cout<<"Tree at 0 "<<TreeNames.at(0)<<std::endl; 
     std::cout<<"Tree at 1 "<<TreeNames.at(1)<<std::endl; 
     if(TreeNames.at(i) == "nominal_Loose")
     {
     TString TreeName  = TString(TreeNames.at(i));
     TChain* mych_data = new TChain (TreeName);
     mych_data->Add(OUTPUTDIR+"/"+path+"/"+SampleName);
     std::cout<<mych_data->GetEntries()<<"  "<<path+"/"+SampleName<<std::endl;
     EventLoop *eventLoop = new EventLoop(mych_data,SampleName,MCDataPeriode,TreeName,WP);
     eventLoop->Loop();
     TDirectory* subdir = outfile->mkdir(TreeNames.at(i).c_str());
     eventLoop->Write(subdir,TreeNames.at(i)); 
     eventLoop->WriteTreesToFile(outfile); 
     delete mych_data;
     delete eventLoop;
     }
  }

  std::cout<<"Finished looping over the events"<<std::endl;  
  outfile->Close();
  return 0;
}

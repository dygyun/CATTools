#ifndef __CINT__

#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>

#endif
// Root
#include "TString.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
//#include "TKey.h"
//#include "TPrint.h"
//#include <exception>
#include <sys/stat.h>
// TopCode
#include <SFIDISOTrigger.h> // SF_ID-ISO-Trigger
#ifndef __CINT__
using namespace std;
void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i inputfile  Input file without .root" << std::endl;
  std::cout << "    -o name in the output file \"h_\"" << std::endl;
  std::cout << "    -s create a file with the systematic uncertainty yields" << std::endl;
  std::cout << "    -tr SF Trigger Uncertainty" << std::endl;
  std::cout << "    -idiso SF ID/ISO Uncertainty" << std::endl;
  std::cout << "    -d Input file directory. Default directory: InputTrees" << std::endl;
  std::cout << "    -h                 displays this help message and exits " << std::endl;
  std::cout << "" << std::endl;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const TString currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d at %X", &tstruct);

  return buf;
}

void print_progress(int TreeEntries, Long64_t ievt)
{
  int step = TreeEntries/50;
  if (ievt%(step) == 0){ 
    float progress=(ievt)/(TreeEntries*1.0);
    int barWidth = 50;
    
    std::cout << "[";
    int pos = barWidth * progress;
    
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }  
}

int main(int argc, const char* argv[]){

  gSystem->Load("libTree");

  bool   _ttbar_cat = false;
  bool   _syst      = false;
  bool   _tr_unc    = false;
  bool   _idiso_unc = false;
  const char * _output   = 0;
  const char * _input    = 0;
  // TopTrees directory
  //const char * _dir      = "/cms/home/dygyun/TopAnal/ForDooyeon/TopCodeljets/ntuple/ABCD/";
  const char * _dir      = "/cms/home/dygyun/TopAnal/ForDooyeon/TopCodeljets/v763/";
  //const char * _dir      = "/Users/dygyun/Top/v746/Javier/";
  const char * _tr       = 0;
  const char * _idiso    = 0;
  const char * _ttbar_id = 0;

  // Arguments used
  //std::set<int> usedargs;
  //Parsing input options
  if(argc == 1){
    display_usage();
    return -1;
  }

  else{
      //Argumet 1 must be a valid input fileName
      for (int i = 1; i < argc; i++){
  if( strcmp(argv[i],"-i") == 0 ){
    _input = argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-d") == 0 ){
    _dir = argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-o") == 0 ){
    _output= argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-s") == 0 ){
    _syst= true;
  }
  if( strcmp(argv[i],"-tr") == 0 ){
    _tr_unc= true;
    _tr= argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-idiso") == 0 ){
    _idiso_unc= true;
    _idiso= argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-cat") == 0 ){
    _ttbar_cat = true;
    _ttbar_id  = argv[i+1];
    i++;
  }
  if( strcmp(argv[i],"-h") == 0 ||
      strcmp(argv[i],"--help") == 0 ){
    display_usage();
    return 0;
  }
      }
  }//else
  if( _input ==0 ){
    std::cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
        << std::endl;
    display_usage();
    return -1;
  }
  
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);
  TString TrUnc(_tr);
  TString IDISOUnc(_idiso);
  TString ttbar_id(_ttbar_id);
  
  TChain theTree("ttbbLepJets/tree"); 
  
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "Signal: ";
  std::cout << fname + ".root" << std::endl;
  std::cout << fdir + fname + ".root" << std::endl;

  theTree.Add(fdir + fname + ".root");
  
  int Event,Run,Channel, GenCat_ID, GoodPV;
  float  GENWeight; 

  float MET,MET_Phi;

  float Lep_px, Lep_py, Lep_pz, Lep_E, Lep_Iso;
  std::vector<float> *Jet_px=0, *Jet_py=0, *Jet_pz=0, *Jet_E=0;
  std::vector<float> *Jet_CSV=0;
  std::vector<float> *PUWeight=0;

  float puWeight_= 0.0;

  /*********************************
           Tree Branches
  **********************************/
  
  theTree.SetBranchAddress( "event",    &Event );
  theTree.SetBranchAddress( "run",      &Run );

  theTree.SetBranchAddress( "PUWeight",  &PUWeight );
  theTree.SetBranchAddress( "channel",   &Channel );

  theTree.SetBranchAddress( "GoodPV",  &GoodPV );

  theTree.SetBranchAddress( "MET",     &MET );
  theTree.SetBranchAddress( "MET_phi", &MET_Phi );

  theTree.SetBranchAddress( "lepton_px", &Lep_px );
  theTree.SetBranchAddress( "lepton_py", &Lep_py );
  theTree.SetBranchAddress( "lepton_pz", &Lep_pz );
  theTree.SetBranchAddress( "lepton_E",  &Lep_E );

  theTree.SetBranchAddress( "jet_px", &Jet_px );
  theTree.SetBranchAddress( "jet_py", &Jet_py );
  theTree.SetBranchAddress( "jet_pz", &Jet_pz );
  theTree.SetBranchAddress( "jet_E",  &Jet_E );

  theTree.SetBranchAddress( "jet_CSV",  &Jet_CSV );

  // Number of Events and Weights (MC@NLO)
  TFile *fileEntries = NULL;
  fileEntries = TFile::Open(fdir + fname + ".root");
  TH1F *h_NumberEvt;
  h_NumberEvt = (TH1F*)fileEntries->Get("ttbbLepJets/EventInfo");

  float NTotal_Event, NTotal_Weight, nNorm_Event;
  NTotal_Event  = h_NumberEvt->GetBinContent(1);
  NTotal_Weight = h_NumberEvt->GetBinContent(2);

  // MCatNLO Weights
  if(fname.Contains("MCatNLO")){
    theTree.SetBranchAddress( "genweight", &GENWeight );
    nNorm_Event = NTotal_Weight;    
  }
  else{
    GENWeight = 1.0;
    nNorm_Event = NTotal_Event;
  }
  /*********************************
             Histograms
  **********************************/
  //Correct Statistical Uncertainty Treatment
  TH1::SetDefaultSumw2(kTRUE);  
  
  TH1F *hPV[4][2];
  TH1F *hMET[4][2], *hMET_Phi[4][2];
  TH2D *hIso_vs_MET[4][2];
  TH1F *hLepPt[4][2], *hLepEta[4][2], *hLepPhi[4][2];
  TH1F *hLepIso[4][2];
  TH1F *hNJets[4][2], *hHT[4][2], *hNBtagJets[4][2];
  TH1F *CSV[4][4][2];

  TH1F *hSFpT[3][2], *hSFpTError[3][2];

  TH1F *hSFIDISO[4][2], *hSFIDISOError[4][2];
  TH1F *hSFTrigger[4][2], *hSFTriggerError[4][2];

  TString namech[2];
  namech[0]="mujets";
  namech[1]="ejets";
  
  TString namecut[4];
  namecut[0]="lepton";
  namecut[1]="4Jets";
  namecut[2]="MET";
  namecut[3]="2btag";
  
  TString titlenamech[2];
  titlenamech[0]="#mu+Jets";
  titlenamech[1]="e+Jets";
  
  for(int j=0; j<4; j++){   // Cut
    for(int i=0; i<2; i++){ // Channel
      hPV[j][i]         = new TH1F("hPV_"+namech[i]+"_"+namecut[j],"PV Distribution  " + titlenamech[i],30,0,30);
      hPV[j][i]->GetXaxis()->SetTitle("PV");      
      hMET[j][i]        = new TH1F("hMET_"+ namech[i]+"_"+namecut[j],"#slash{E}_{T} " + titlenamech[i],20,0,200);
      hMET[j][i]->GetXaxis()->SetTitle("#slash{E}_{T}[GeV]");      
      hMET_Phi[j][i]    = new TH1F("hMET_Phi_"+ namech[i]+"_"+namecut[j],"#Phi_{#slash{E}_{T}} " + titlenamech[i],160,-4,4);
      hMET_Phi[j][i]->GetXaxis()->SetTitle("#Phi_{#slash{E}_{T}}[rad]");      
      hIso_vs_MET[j][i]    = new TH2D("hIso_vs_MET"+ namech[i]+"_"+namecut[j],"Iso_{#slash{E}_{T}} " + titlenamech[i],100,0,300,100,0,4.);
      hIso_vs_MET[j][i]->GetXaxis()->SetTitle("PF Isolation vs {#slash{E}_{T}}");      
      
      hLepPt [j][i]    = new TH1F("hLepPt_"  + namech[i] + "_" + namecut[j], "Lepton p_{T} " + titlenamech[i],20,0.0,200.0);
      hLepPt[j][i]->GetXaxis()->SetTitle("Lepton p_{T}[GeV]");      
      hLepEta[j][i]    = new TH1F("hLepEta_" + namech[i] + "_" + namecut[j], "#eta_{Lep} " + titlenamech[i],12,0.0,2.2);
      hLepEta[j][i]->GetXaxis()->SetTitle("Lepton Abs(#eta)");      
      hLepPhi[j][i]    = new TH1F("hLepPhi_" + namech[i] + "_" + namecut[j], "#phi_{Lep} " + titlenamech[i],16,0.0,3.2);
      hLepPhi[j][i]->GetXaxis()->SetTitle("lepton #Phi[rad]");      
      hLepIso[j][i]    = new TH1F("hLepIso"  + namech[i] + "_" + namecut[j], "Lepton PF Isolation" + titlenamech[i],20,0.0,2.0);
      hLepIso[j][i]->GetXaxis()->SetTitle("Lepton Isolation");      
      
      hNJets[j][i]      = new TH1F("hNJets_" + namech[i] + "_" + namecut[j], "Jet multiplicity " + titlenamech[i],9,-0.5,8.5);
      hNJets[j][i]->GetXaxis()->SetTitle("Number of jets");      
      hNJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNJets[j][i]->GetXaxis()->SetBinLabel(4,"3");
      hNJets[j][i]->GetXaxis()->SetBinLabel(5,"4");
      hNJets[j][i]->GetXaxis()->SetBinLabel(6,"5");
      hNJets[j][i]->GetXaxis()->SetBinLabel(7,"6");
      hNJets[j][i]->GetXaxis()->SetBinLabel(8,"7");
      hNJets[j][i]->GetXaxis()->SetBinLabel(9,"#geq 8");
      
      hNBtagJets[j][i]  = new TH1F("hNBtagJets_"+namech[i]+"_"+namecut[j],"b-tag jet multiplicity " + titlenamech[i],4,-0.5,3.5);
      hNBtagJets[j][i]->GetXaxis()->SetTitle("Number of b-jets");      
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(4,"#geq 3");
      
      hHT[j][i]         = new TH1F("hHT_"+namech[i]+"_"+namecut[j],"H_{T} " + titlenamech[i],100,0,600);
      hHT[j][i]->GetXaxis()->SetTitle("HT[GeV]");      

      /***************************
          SF(ID,ISO & Trigger)
      ***************************/
      hSFIDISO[j][i]           = new TH1F("hSFIDISO_"+namech[i]+"_"+namecut[j],"SF_{ID,ISO} " + titlenamech[i],400,0.8,1.2);    
      hSFIDISOError[j][i]      = new TH1F("hSFIDISOError_"+namech[i]+"_"+namecut[j],"#Delta SF_{ID,ISO} " + titlenamech[i],400,0,0.05); 
      hSFTrigger[j][i]         = new TH1F("hSFTrigger_"+namech[i]+"_"+namecut[j],"SF^{Trigger} " + titlenamech[i],400,0.8,1.2);    
      hSFTriggerError[j][i]    = new TH1F("hSFTriggerError_"+namech[i]+"_"+namecut[j],"#Delta SF^{Trigger} " + titlenamech[i],400,0,0.05);
            
      /***************************
            SF(pT Reweight)
      ***************************/
      hSFpT[j][i]           = new TH1F("hSFpT_"+namech[i]+"_"+namecut[j],"SF_{pT} " + titlenamech[i],500,0.4,1.4);    
      hSFpTError[j][i]      = new TH1F("hSFpTError_"+namech[i]+"_"+namecut[j],"#Delta SF_{pT} " + titlenamech[i],400,0,0.05); 
      
      
      TString jetn[4];
      jetn[0]= "Jet-0"; 
      jetn[1]= "Jet-1"; 
      jetn[2]= "Jet-2"; 
      jetn[3]= "Jet-3"; 
      
      for(int ij=0; ij<4; ij++){
  CSV[ij][j][i]         = new TH1F("hCSV_" + jetn[ij] + "_" + namech[i] + "_" + namecut[j],"CSV " + jetn[ij] + " " + titlenamech[i],20,0,1);
      CSV[ij][j][i]->GetXaxis()->SetTitle("CSVv2");      
      }

    }//for(i)
  }//for(j)
  

  TStopwatch sw;
  sw.Start(kTRUE);

  ///////////////////////////////////////
  // Please, IGNORE. Temporal solution //
  ///////////////////////////////////////
  TCanvas *mydummycanvas=new TCanvas();// 
  ///////////////////////////////////////
    
  /************************
     SF Parametrization
  *************************/

 // TString fSFdir = "/cms/home/dygyun/TopAnal/ForDooyeon/TopCodeljets/ntuple/ScaleFactors/";
  TString fSFdir = "/cms/home/brochero/ScaleFactors/";
  
  TH2F *hmuIDISOSF, *hmuTriggerSF;
  TH2F *heIDISOSF,  *heTriggerSF;

  // Lepton SFs: ID and ISO with stat. + syst. Errors
  TFile *MuSF = TFile::Open(fSFdir + "SF_muon_IDISO_13TeV_v2.root"); 
  TFile *ElSF = TFile::Open(fSFdir + "SF_electron_IDISO_13TeV_v2.root"); 

  if(!MuSF || !ElSF){
    std::cerr << "ERROR [SF]: Could not open SF files!!!"  << std::endl;
    std::exit(0);
  }

  hmuIDISOSF = (TH2F*) MuSF->Get("GlobalSF")->Clone("muIDISOSF");
  hmuTriggerSF = (TH2F*) MuSF->Get("TriggerSF")->Clone("muTriggerSF"); 
  if(!hmuIDISOSF || !hmuTriggerSF){
    std::cerr << "ERROR [MuonSF]: Could not find histogram for SF reweighting" << std::endl;
  }

  heIDISOSF = (TH2F*) ElSF->Get("GlobalSF")->Clone("eIDISOSF");
  heTriggerSF = (TH2F*) ElSF->Get("TriggerSF")->Clone("eTriggerSF"); 
  if(!heIDISOSF || !heTriggerSF){
    std::cerr << "ERROR [ElectronSF]: Could not find histogram for SF reweighting" << std::endl;
  }

  // Trigger and ID-ISO uncertainties
  if(_idiso_unc){
    if(IDISOUnc=="Up")        fname += "_SYS_IDISO_Up";      
    else if(IDISOUnc=="Down") fname += "_SYS_IDISO_Down";
    else if(IDISOUnc=="Nom")  fname += "_SYS_IDISO_Nom";
  } // if(_idiso_unc)
  
  if(_tr_unc){
    if(TrUnc=="Up")        fname += "_SYS_Trigger_Up";
    else if(TrUnc=="Down") fname += "_SYS_Trigger_Down";
    else if(TrUnc=="Nom")  fname += "_SYS_Trigger_Nom";
  }// if(_tr_unc) 


  // Number de events for <pT Reweight>
  //          [Cut][Channel]
  float SF_pTweight[4][3]={0,0,0,0,
         0,0,0,0,
         0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  int AccEvent[4][3]={0,0,0,0,
          0,0,0,0,
          0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  float EffEvent[4][3]={0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0};
  
  /************************
    Normalization Weights
  *************************/
  
  float Lumi=2260.; //pb-1 From the JSON file Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt
  float NormWeight = 0.0;
  // NormWeight = Lumi*(1.0/N_Gen_events)*(Xsec)*(Br)
  if(fname.Contains("DYJets"))         NormWeight = Lumi * (1.0/nNorm_Event) * (6025.2);  // [pb]
  if(fname.Contains("DYJets_10to50"))  NormWeight = Lumi * (1.0/nNorm_Event) * (18610);// [pb]
  if(fname.Contains("WJets"))          NormWeight = Lumi * (1.0/nNorm_Event) * (61526.7); // [pb]
  if(fname.Contains("SingleTop_tW"))   NormWeight = Lumi * (1.0/nNorm_Event) * (35.6);    // [pb]
  if(fname.Contains("SingleTbar_tW"))  NormWeight = Lumi * (1.0/nNorm_Event) * (35.6);    // [pb]
  if(fname.Contains("SingleTop_t"))    NormWeight = Lumi * (1.0/nNorm_Event) * (44.33);  // [pb]
  if(fname.Contains("SingleTbar_t"))   NormWeight = Lumi * (1.0/nNorm_Event) * (26.38);   // [pb]
  if(fname.Contains("WW"))             NormWeight = Lumi * (1.0/nNorm_Event) * (118.7);   // [pb]
  if(fname.Contains("WZ"))             NormWeight = Lumi * (1.0/nNorm_Event) * (47.13);    // [pb]
  if(fname.Contains("ZZ"))             NormWeight = Lumi * (1.0/nNorm_Event) * (16.523);    // [pb]

  if(fname.Contains("TT_powheg"))      NormWeight = Lumi * (1.0/nNorm_Event) * (831.76); // [pb] Br = (leptonic) * Hadronic = (0.1086*3) * (0.67)
  if(fname.Contains("TTJets_aMC"))     NormWeight = Lumi * (1.0/nNorm_Event) * (831.76); // * (0.1086*3.0*3.0); // Br correction  
  if(fname.Contains("TTJets_MG5"))     NormWeight = Lumi * (1.0/nNorm_Event) * (831.76); //

  if(fname.Contains("Run2015"))       NormWeight = 1.0;  
  if(fname.Contains("QCD_Pt-20to30_MuEnriched"))        NormWeight = Lumi * (1.0/nNorm_Event) * (2960198.4) * (0.0053);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-30to50_MuEnriched"))        NormWeight = Lumi * (1.0/nNorm_Event) * (1652471.46) * (0.0053);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-50to80_MuEnriched"))        NormWeight = Lumi * (1.0/nNorm_Event) * (437504.1) * (0.02276);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-80to120_MuEnriched"))       NormWeight = Lumi * (1.0/nNorm_Event) * (106033.6648) * (0.03844);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-120to170_MuEnriched"))      NormWeight = Lumi * (1.0/nNorm_Event) * (106033.6648) * (0.05362);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-170to300_MuEnriched"))      NormWeight = Lumi * (1.0/nNorm_Event) * (8654.49315) * (0.07335);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-300to470_MuEnriched"))      NormWeight = Lumi * (1.0/nNorm_Event) * (797.35269) * (0.10196);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-470to600_MuEnriched"))      NormWeight = Lumi * (1.0/nNorm_Event) * (79.02553776) * (0.12242);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-600to800_MuEnriched"))      NormWeight = Lumi * (1.0/nNorm_Event) * (25.09505908) * (0.12242);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-800to1000_MuEnriched"))     NormWeight = Lumi * (1.0/nNorm_Event) * (4.707368272) * (0.14552);   // [pb] (cross section) * (Filter Eff)
  if(fname.Contains("QCD_Pt-1000toInf_MuEnriched"))     NormWeight = Lumi * (1.0/nNorm_Event) * (1.62131692) * (0.15544);   // [pb] (cross section) * (Filter Eff)


  std::cout << "-----------------------                                 -------------------------" << std::endl;
  std::cout << "Number of Events     = " << nNorm_Event << std::endl;
  std::cout << "Normalization Factor = " << NormWeight  << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;

  /********************************
      pT Reweight; <SF_pT> 
   *******************************/
  // <SF_pT> per channel at single lepton level 
  float SF_tPt_mean[4][3];
  SF_tPt_mean[0][0]= 1.0;
  SF_tPt_mean[0][1]= 1.0;
  SF_tPt_mean[0][2]= 1.0;
  SF_tPt_mean[1][0]= 1.0;
  SF_tPt_mean[1][1]= 1.0;
  SF_tPt_mean[1][2]= 1.0;
  SF_tPt_mean[2][0]= 1.0;
  SF_tPt_mean[2][1]= 1.0;
  SF_tPt_mean[2][2]= 1.0;
  SF_tPt_mean[3][0]= 1.0;
  SF_tPt_mean[3][1]= 1.0;
  SF_tPt_mean[3][2]= 1.0;
  
  float SF_pT_TotalMean; // Over all cuts and channels
  float SF_n = 0.0;
  for(int SF_cut= 0; SF_cut<2; SF_cut++){
    for(int SF_ch= 0; SF_ch<3; SF_ch++){
      SF_pT_TotalMean +=  SF_tPt_mean[SF_cut][SF_ch];
      SF_n ++;
    }
  }
  if(SF_n != 0.0) SF_pT_TotalMean = SF_pT_TotalMean/SF_n; 
  
  /********************************
             Event Loop
  ********************************/
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;

  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
  //  cout << "--- test -- " << endl;
    theTree.GetEntry(ievt);  
    print_progress(theTree.GetEntries(), ievt);


    // MCatNLO GEN Weights (For MC@NLO)
    //PUWeight = PUWeight * GENWeight; // original //
    puWeight_ = PUWeight->at(0) * GENWeight; // //
    // Normalization Weight
    puWeight_ = PUWeight->at(0) * NormWeight;
    
    int NJets,NBtagJets;
    
    TLorentzVector Lep;
    std::vector<TLorentzVector> Jet;
    std::vector<TLorentzVector> bJet;
         
    Lep.SetPxPyPzE(Lep_px,Lep_py,Lep_pz,Lep_E);
    if(Lep.Pt() < 30)  continue; // Lep pT >30GeV

    // Jets 
    NJets     = 0;
    NBtagJets = 0;

    for(int ijet=0; ijet < Jet_px->size(); ijet++){
      
      TLorentzVector jet;
      jet.SetPxPyPzE((*Jet_px)[ijet],(*Jet_py)[ijet],(*Jet_pz)[ijet],(*Jet_E)[ijet]);
      
      if(jet.Pt()>25){ // Jet pT > 25GeV
	
	Jet.push_back(jet);
	NJets++; // Number of jets
	
	// b-tagging WP from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging
	if((*Jet_CSV)[ijet] > 0.800){ // CSVM. Luca b-tagging code to apply SF?
	  bJet.push_back(jet);
	  NBtagJets++; // Number of b-tagged jets
	} // if(b-tag)
      } // if(Jet_pT)
    }// for(jets)
    
    
  /*******************************************
   Trigger,ID & ISO Scale Factors/bin(Pt,Eta)
  *******************************************/    
    std::vector<float> SF_ID_ISO_Tr;
    
    if (fname.Contains("Data")){
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO_Tr    [0] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO       [1] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO_Error [2] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Tr           [3] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Tr_Error     [4] 
    }
    
    else {
      SFIDISOTrigger(SF_ID_ISO_Tr,
		     Lep, Channel,
		     hmuIDISOSF, hmuTriggerSF,
		     heIDISOSF,  heTriggerSF);
      
      if(_idiso_unc){
	if     (IDISOUnc == "Up")   puWeight_ = puWeight_ * (SF_ID_ISO_Tr[1] + SF_ID_ISO_Tr[2]);	
	else if(IDISOUnc == "Down") puWeight_ = puWeight_ * (SF_ID_ISO_Tr[1] - SF_ID_ISO_Tr[2]);	
	else if(IDISOUnc == "Nom")  puWeight_ = puWeight_ * (SF_ID_ISO_Tr[1]); // 
      } // if(_idiso_unc)
      
      else if(_tr_unc){
	if     (TrUnc=="Up")   puWeight_=puWeight_*(SF_ID_ISO_Tr[3] + SF_ID_ISO_Tr[4]);			
	else if(TrUnc=="Down") puWeight_=puWeight_*(SF_ID_ISO_Tr[3] - SF_ID_ISO_Tr[4]);		
	else if(TrUnc=="Nom")  puWeight_=puWeight_*(SF_ID_ISO_Tr[3]);	// 
      }// if(_tr_unc) 
      
      // If the electron is in the transition region, SF = 1
      else  if(SF_ID_ISO_Tr[0] != 0.0) puWeight_=puWeight_*(SF_ID_ISO_Tr[0]); //
      
    }// else(Contain("Data"))
    

    /***************************
            Selection
    ***************************/

    int                                    cut = 0; // Single Lepton (from Tree)
    if(NJets>3)                            cut = 1; // + 4Jets 
    if(NJets>3 && MET>30.0)                cut = 2; // + MET
    if(NJets>3 && MET>30.0 && NBtagJets>1) cut = 3; // + 2 Btag

    /***************************
        ttbar Categorization
     ***************************/
//    if (_ttbar_cat && !ttbar_category(ttbar_id, GenCat_ID))  cut = -1;

    /***************************
          Loop over cuts
    ***************************/
    for(int icut = 0; icut < (cut+1); icut++){
            
      /************************************************************
        pT reweight (Only for ttbar signal)
	It shouldn't be applied to the central value. 
        It is divided over the <SF_pT> (normalization)  
      *************************************************************/
      if(fname.Contains("ttbar")){	
	// // pT Top Reweight
	// TLorentzVector t,tbar;
	// t.SetPxPyPzE(tPx,tPy,tPz,tE);
	// tbar.SetPxPyPzE(tbarPx,tbarPy,tbarPz,tbarE);
	
	// // pT Reweight: Only for Systematic Uncertainty
	// // Definition from pT reweight Twiki
	// // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting and 
	// // https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematicsRun1#pt_top_Reweighting
	// float SF_tPt=sqrt( exp(0.222-0.00197*TMath::Min(t.Pt(), 400.)) * exp(0.222-0.00197*TMath::Min(tbar.Pt(), 400.)) );
	// hSFpT[icut][Channel]->Fill(SF_tPt,PUWeight); 	
	
	// if(fname.Contains("TTJets_Dilep_SYS_pTReweight_Up") || fname.Contains("TTJets_Dilep_SYS_pTReweight_Down")) PUWeight=PUWeight*(SF_tPt/SF_pT_TotalMean); 
	// SF_pTweight[icut][Channel]+=SF_tPt;
	SF_pTweight[icut][Channel] += 1.0;	
      } // if(TTbar) 
      /*************************************************************/
      
      /*******************
        Fill Histograms
      ******************/

      hSFIDISO[icut][Channel]->Fill(SF_ID_ISO_Tr[1],puWeight_);
      hSFIDISOError[icut][Channel]->Fill(SF_ID_ISO_Tr[2],puWeight_);
      hSFTrigger[icut][Channel]->Fill(SF_ID_ISO_Tr[3],puWeight_);
      hSFTriggerError[icut][Channel]->Fill(SF_ID_ISO_Tr[4],puWeight_);
    
      /******************
          Acceptace
      ******************/
      AccEvent[icut][Channel]++;
      EffEvent[icut][Channel]= EffEvent[icut][Channel] + puWeight_;


      hPV[icut][Channel]->Fill(GoodPV,puWeight_);

      /******************
        Kinematic Var.
      ******************/
      
      hMET[icut][Channel]->Fill(MET,puWeight_);
      hMET_Phi[icut][Channel]->Fill(fabs(MET_Phi),puWeight_);
      
      hLepPt[icut][Channel]->Fill(Lep.Pt(),puWeight_);
      hLepEta[icut][Channel]->Fill(fabs(Lep.Eta()),puWeight_);
      hLepPhi[icut][Channel]->Fill(fabs(Lep.Phi()),puWeight_);
      
      /******************
          Jets Var.
      ******************/

      hNJets[icut][Channel]->Fill(NJets,puWeight_); 
      hNBtagJets[icut][Channel]->Fill(NBtagJets,puWeight_);

      for(int ijet=0; ijet < Jet_px->size(); ijet++){
	if (ijet < 4) CSV[ijet][icut][Channel]->Fill((*Jet_CSV)[ijet],puWeight_);
      }
    }//for(icuts)     
    
    Jet.clear();
    bJet.clear();

  }//for(events)
  

  /***************************
      TTbar pT Reweight
  ***************************/

  if(fname.Contains("ttbar")){
    
    for(int icut=0; icut<4;icut++){
      for(int ich=0; ich<2;ich++){
  
  // SF_pT_reweight mean value 
  SF_pTweight[icut][ich]=SF_pTweight[icut][ich]/AccEvent[icut][ich]; 
  std::cout << "<SF_pT[cut=" << icut << "][ch=" <<ich << "]>= " << SF_pTweight[icut][ich] << std::endl; 
      } // for(icut)
    } // for(ich)
  } // if(TT)


  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  

  //Acceptance-Efficiency
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-mu+Jets events:" << std::endl;
  std::cout << "lepton: "  << AccEvent[0][0] << std::endl;
  std::cout << "4 Jets: "  << AccEvent[1][0] << std::endl;
  std::cout << "MET:    "  << AccEvent[2][0] << std::endl;
  std::cout << "2btag:  "  << AccEvent[3][0] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-mu+Jets events:" << std::endl;
  std::cout << "lepton: "  << EffEvent[0][0] << " +/- " << sqrt(AccEvent[0][0])*EffEvent[0][0]/AccEvent[0][0] << std::endl;
  std::cout << "4 Jets: "  << EffEvent[1][0] << " +/- " << sqrt(AccEvent[1][0])*EffEvent[1][0]/AccEvent[1][0] << std::endl;
  std::cout << "MET:    "  << EffEvent[2][0] << " +/- " << sqrt(AccEvent[2][0])*EffEvent[2][0]/AccEvent[2][0] << std::endl;
  std::cout << "2btag:  "  << EffEvent[3][0] << " +/- " << sqrt(AccEvent[3][0])*EffEvent[3][0]/AccEvent[3][0] << std::endl;
  
  std::cout << "--------  Acceptace  --------" << std::endl;
  std::cout << "Number of RAW-e+Jets events:" << std::endl;
  std::cout << "lepton: "   << AccEvent[0][1] << std::endl;
  std::cout << "4 Jets: "   << AccEvent[1][1] << std::endl;
  std::cout << "MET:    "   << AccEvent[2][1] << std::endl;
  std::cout << "2 btag: "   << AccEvent[3][1] << std::endl;

  std::cout << "--------  Efficiency  --------" << std::endl;
  std::cout << "Number of Weigthed-e+Jets events: " << std::endl;
  std::cout << "lepton: "   << EffEvent[0][1] << " +/- " << sqrt(AccEvent[0][1])*EffEvent[0][1]/AccEvent[0][1] << std::endl;
  std::cout << "4 Jets: "   << EffEvent[1][1] << " +/- " << sqrt(AccEvent[1][1])*EffEvent[1][1]/AccEvent[1][1] << std::endl;
  std::cout << "MET:    "   << EffEvent[2][1] << " +/- " << sqrt(AccEvent[2][1])*EffEvent[2][1]/AccEvent[2][1] << std::endl;
  std::cout << "2btag:  "   << EffEvent[3][1] << " +/- " << sqrt(AccEvent[3][1])*EffEvent[3][1]/AccEvent[3][1] << std::endl;
  //Output Dir
  TString dirname="TopResults_wmass";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);
  
  // Sample name identification
  TString samplename="";
  bool matchsamplename=false;
  
  for(int i=0; i<fname.Sizeof(); i++){
    if (i>2){
      if (fname[i-3]=='-' && 
    fname[i-2]=='1' && 
    fname[i-1]=='_') matchsamplename=true;
    }
    if (matchsamplename) samplename.Append(fname[i]);
  }
  
  if(!_syst){
    // Yields
    // make a dir if it does not exist!!
    TString diryieldsname = dirname + "/Yields_" + hname;
    struct stat st;
    if(stat(diryieldsname,&st) != 0) system("mkdir " + diryieldsname);
    
    TString Yieldfile = diryieldsname + "/" + samplename.Data() + ".h";
    //Yieldfile += ".h";
    // Option a = append: Open file for output at the end of a file.
    // Option w = write: Create an empty file for output operations. 
    FILE* fyields = fopen(Yieldfile, "w"); 

    fprintf(fyields,"\n///////////////////////////////////////////////////////////////////////////////// \n\n");
    fprintf(fyields,"// %s Sample on %s \n", (fname + ".root").Data() , currentDateTime().Data());
    fprintf(fyields,"// %s version \n", hname.Data());
    fprintf(fyields," float  %s[16][3][2][4];//[systematic][variation][channel][Cut] \n",      samplename.Data());
    fprintf(fyields," float  err_%s[16][3][3][4]; //[systematic][variation][channel][Cut] \n", samplename.Data());
    fprintf(fyields,"// Systematic: [0]=Trigger [1]=ID-ISO [2]=LES \n");
    fprintf(fyields,"// Systematic: [3]=JES [4]=JER [5]=b-tag [6]=PileUp \n");
    fprintf(fyields,"// Systematic: [7]=Scale [8]=Matching [9]=pTreweight [9]=PDF \n");
    fprintf(fyields,"// Systematic: [10]=DY-DD [11]=Non-W/Z [15]=Nom \n");
    fprintf(fyields,"// Variation: [0]=Up [1]=Down [2]=Nom \n");
    fprintf(fyields,"// Channel: [0]=mumu [1]=ee [2]=mue \n");
    fprintf(fyields,"// Cut: [0]=Dielpton [1]=Jets+Zveto [2]=MET [3]=btag \n");

    for(int ch=0;ch<2;ch++){
      for(int cut=0;cut<4;cut++){
      fprintf(fyields,"%s[15][2][%i][%i] = %.3f ; \n", samplename.Data(), ch, cut, EffEvent[cut][ch]);
      if(AccEvent[cut][ch]!=0.0) fprintf(fyields,"err_%s[15][2][%i][%i] = %.3f ; \n", samplename.Data(), ch, cut, sqrt(AccEvent[cut][ch])*EffEvent[cut][ch]/AccEvent[cut][ch]);
      else fprintf(fyields,"err_%s[15][2][%i][%i] = 0.0 ; \n", samplename.Data(), ch, cut);
      }
    }
    fclose(fyields);
    
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Yields saved into " << Yieldfile << " file" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;    
  } 
  
  if (_syst){    

    // Name change
    TString processname = samplename;
    TString systname    = "";
    bool matchsystname  = false; 
    
    for(int i=0;i<samplename.Sizeof();i++){
      if (i>3){
  if (samplename[i-4]=='S' && 
      samplename[i-3]=='Y' && 
      samplename[i-2]=='S' && 
      samplename[i-1]=='_'){
    matchsystname = true;
  }
      }
      if (matchsystname) systname.Append(samplename[i]);
    }
    
    processname.ReplaceAll(("_SYS_" + systname).Data(), "");  
    
    int variation  = -999;
    int systsource = -999;
    
    if     (systname.Contains("Up"))   variation = 0;
    else if(systname.Contains("Down")) variation = 1;
    else if(systname.Contains("Nom"))  variation = 2;

    if(systname.Contains("Trigger"))         systsource = 0;
    else if(systname.Contains("IDISO"))      systsource = 1;
    else if(systname.Contains("LES"))        systsource = 2;
    else if(systname.Contains("JES"))        systsource = 3;
    else if(systname.Contains("JER"))        systsource = 4;
    else if(systname.Contains("btag"))       systsource = 5;
    else if(systname.Contains("PU"))         systsource = 6;
    else if(systname.Contains("Scale"))      systsource = 7;
    else if(systname.Contains("Matching"))   systsource = 8;
    else if(systname.Contains("pTReweight")) systsource = 9;
    else if(systname.Contains("Powheg"))     systsource = 10;

    // Systematic Uncertainty Estimations
    // make a dir if it does not exist!!
    TString dirSysyieldsname = dirname + "/SysYields_" + hname;
    struct stat st;
    if(stat(dirSysyieldsname,&st) != 0) system("mkdir " + dirSysyieldsname);
    
    TString Syst_Yieldfile = dirSysyieldsname + "/" + samplename.Data() + ".h";
    FILE* fSys = fopen(Syst_Yieldfile, "w");        
    
    fprintf(fSys,"\n///////////////////////////////////////////////////////////////////////////////// \n\n");
    fprintf(fSys,"// %s Sample on %s \n", (fname + ".root").Data(), currentDateTime().Data());
    fprintf(fSys,"// %s version \n", hname.Data());
    fprintf(fSys,"// float  %s[16][3][3][4];//[systematic][variation][channel][Cut] \n", processname.Data());
    fprintf(fSys,"// float  err_%s[16][3][3][4]; //[systematic][variation][channel][Cut] \n", processname.Data());
    fprintf(fSys,"// Systematic: [0]=Trigger [1]=ID-ISO [2]=LES \n");
    fprintf(fSys,"// Systematic: [3]=JES [4]=JER [5]=b-tag [6]=PileUp \n");
    fprintf(fSys,"// Systematic: [7]=Scale [8]=Matching [9]=pTreweight [9]=PDF \n");
    fprintf(fSys,"// Systematic: [10]=DY-DD [11]=Non-W/Z [15]=Nom \n");
    fprintf(fSys,"// Variation: [0]=Up [1]=Down [2]=Nom \n");
    fprintf(fSys,"// Channel: [0]=mumu [1]=ee [2]=mue \n");
    fprintf(fSys,"// Cut: [0]=Dielpton [1]=Jets+Zveto [2]=MET [3]=btag \n");

    for(int ch=0;ch<3;ch++){
      for(int cut=0;cut<4;cut++){
  fprintf(fSys,"%s     [%i][%i][%i][%i] = %.3f ; \n", processname.Data(), systsource, variation, ch, cut, EffEvent[cut][ch]);
  if(AccEvent[cut][ch]!=0.0) fprintf(fSys,"err_%s [%i][%i][%i][%i] = %.3f ; \n", processname.Data(), systsource, variation, ch, cut, sqrt(AccEvent[cut][ch])*EffEvent[cut][ch]/AccEvent[cut][ch]);
  else fprintf(fSys,"err_%s [%i][%i][%i][%i] = 0.0 ; \n", processname.Data(), systsource, variation, ch, cut);
      }
    }
    fclose(fSys);
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Yields saved for Syst. estimation into " << Syst_Yieldfile << " file" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;    
  }
  // --- Write histograms
  TString outfname=dirname + "/hSF-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  

  for(int j=0; j<4; j++){
    for(int i=0; i<2; i++){
      
      hPV[j][i]->Write();
      
      hMET[j][i]->Write();
      hMET_Phi[j][i]->Write();
      hIso_vs_MET[j][i]->Write();      

      hLepPt[j][i]->Write();
      hLepEta[j][i]->Write();
      hLepPhi[j][i]->Write();
      
      hNJets[j][i]->Write();
      hNBtagJets[j][i]->Write();            
      
      for(int ij=0; ij<4; ij++){
        CSV[ij][j][i]->Write();
      }
      
      hSFpT[j][i]->Write();
      hSFpTError[j][i]->Write();
      
      hSFIDISO[j][i]->Write();
      hSFIDISOError[j][i]->Write();
      hSFTrigger[j][i]->Write();
      hSFTriggerError[j][i]->Write();
    }//for(i)
  }//for(j)
  std::cout << "File saved as " << outfname << std::endl;
}

#endif

//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <csignal>
#include <numeric>
#include <climits>

/////////////////////////   ROOT   //////////////////////////
#include <TROOT.h>

/////////////////////////   RAT    //////////////////////////
#include <RAT/DS/EV.hh>

/////////////////////////   USER   //////////////////////////
#include "utils.hh"
#include "FillBackRATMC.hh"

#include "Analyzer.hh"
#include "AnalyzerFunctions.hh"

#include "ProgressBar.hpp"

using namespace std;

int main(int argc, char *argv[]) {

  // Get Signal if user wants to interrupt loop
  EoF=0;
  signal(SIGINT,Interrupt);


  // Create TApp
  TApplication theApp("App", &argc, argv);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####            PARSE AND SET INPUT PARAMS             #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  // Input parameters
  // There's default value chosen by ternaries below
  string inputRATName;
  string inputFLATName;

  auto User_nEvts = INT_MIN;
  auto User_iEvt = INT_MIN;

  auto isVerbose = false;

  ProcessArgs(&theApp, &isVerbose,
			  &User_nEvts, &User_iEvt,
			  &inputRATName, &inputFLATName);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE ANALYZER                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto FileAnalyzer = new Analyzer(inputRATName.c_str());
  unsigned long nEvts = User_nEvts > INT_MIN && User_nEvts < FileAnalyzer->GetNEvts() ?
						User_nEvts : FileAnalyzer->GetNEvts();
  unsigned long iEvt = SetDefValue(User_iEvt, 0);
  nEvts = iEvt > 0 ? iEvt + nEvts : nEvts;
  nEvts = nEvts > FileAnalyzer->GetNEvts() ? FileAnalyzer->GetNEvts() : nEvts;
  ProgressBar progressBar(nEvts, 70);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####           LOOP AND FILL MAP USING MCID            #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  map<int, int> mEntry;

  if(isVerbose)
    cout << "Recovering EV object for each events" << endl;

  for(iEvt; iEvt<nEvts; iEvt++) {

	// record the tick
	++progressBar;

	auto mc = GetRATMCOnEvt(FileAnalyzer, iEvt);
	mEntry[mc->GetID()] = iEvt;

    if(isVerbose)
      progressBar.display();

  }

  progressBar.done();

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####      FILL RECON DATA TO ORIGINAL RAT MC           #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  FileAnalyzer->GetTreeMc()->SetTreeIndex(0);

  auto tOutput = new TTree("OffT", "Offline Recon EV Tree");
  // auto PF = new RAT::DS::PathFit();
  // auto bPathFit = tOutput->Branch("OffEV", "RAT::DS::PathFit", &PF);
  auto OffEV = new RAT::DS::EV();
  auto bPathFit = tOutput->Branch("OffEV", "RAT::DS::EV", &OffEV);

  auto tFile = TFile::Open(inputFLATName.c_str());

  double X, Y, Z, T, Theta, Phi;
  Long64_t MCID;
  auto TRecon = SetFlatTreeReader(tFile,
								  X, Y, Z,
								  T,
								  Theta, Phi,
								  MCID);

  auto nReconEntries = TRecon->GetEntries();

  nReconEntries = nReconEntries >= nEvts ? nEvts : nReconEntries;

  ProgressBar progressBarRecon(nReconEntries, 70);

  if(isVerbose)
    cout << "Filling EV objects with Reconstructed variables" << endl;

  for(auto iEntry=0; iEntry<nReconEntries; iEntry++){

	// record the tick
	++progressBarRecon;

	tFile->cd();
	TRecon->GetEntry(iEntry);

	FileAnalyzer->GetFmc()->cd();
	OffEV = GetRATEVOnEvt(FileAnalyzer, mEntry[MCID]);

	OffEV->GetPathFit()->SetTime0(T);
	OffEV->GetPathFit()->SetTime(T);

	OffEV->GetPathFit()->SetPos0(TVector3(X, Y, Z));

	TVector3 dir(0.,0.,1.);
	dir.SetMagThetaPhi(1, Theta, Phi);
	OffEV->GetPathFit()->SetDirection(dir);

	bPathFit->Fill();
	// OffEV->Clear();

	if(isVerbose)
	  progressBarRecon.display();

  }

  progressBarRecon.done();

  delete TRecon;
  tFile->Close();

  auto OffEVFile = new TFile("OffEVFile.root", "RECREATE");
  OffEVFile->cd();
  tOutput->Write();
  OffEVFile->Close();

  FileAnalyzer->GetTreeMc()->AddFriend(tOutput);
  FileAnalyzer->GetFmc()->Write();
  FileAnalyzer->GetFmc()->Close();

  /////////////////////////
  // ...

  if(gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}

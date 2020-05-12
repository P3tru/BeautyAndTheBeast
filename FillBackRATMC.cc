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
#include <RAT/DS/Root.hh>

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

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####           LOOP AND FILL MAP USING MCID            #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  map<int, int> mMCID
	  = GetMAssociatedIDAndEvt(FileAnalyzer, iEvt, nEvts, isVerbose);
  delete FileAnalyzer;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####           LOOP AND FILL MAP USING MCID            #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  std::vector< pair<int, FlatRecon> > vpFlat
	  = GetVAssociatedIDAndRecon(inputFLATName.c_str(), isVerbose);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####      FILL RECON DATA TO ORIGINAL RAT MC           #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto tOutput = new TTree("OffT", "T With Offline Recon");
  auto OffDS = new RAT::DS::Root();
  tOutput->Branch("OffDS", "RAT::DS::Root", &OffDS);

  if(isVerbose)
	cout << "Writing recon info into DS object" << endl;

  ProgressBar progressBarRecon(nEvts, 70);

  for(iEvt; iEvt<nEvts; iEvt++){

	// record the tick
	++progressBarRecon;

	auto fMC = TFile::Open(inputRATName.c_str());
	auto fMCTree = dynamic_cast<TTree *>(fMC->Get("T"));
	auto BufDS = new RAT::DS::Root();
	fMCTree->SetBranchAddress("ds", &BufDS);
	fMCTree->GetEntry(iEvt);

	fMC->Close();
	delete fMC;

	AddReconInfo(BufDS,
				 0,
				 find_if(vpFlat.begin(), vpFlat.end(),
						 [&](const pair<int, FlatRecon>& p){return p.first == mMCID[iEvt];})->second);

	tOutput->SetBranchAddress("OffDS", &BufDS, nullptr);
	int errcode = tOutput->Fill();
	tOutput->SetBranchAddress("OffDS", &OffDS, nullptr);


	if(isVerbose)
	  progressBarRecon.display();

  }

  if(isVerbose)
	progressBarRecon.done();

  // delete fA;

  auto OffEVFile = new TFile("OffEVFile.root", "RECREATE");
  tOutput->Write();
  OffEVFile->Close();

  /////////////////////////
  // ...

  if(gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}

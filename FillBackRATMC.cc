//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <csignal>
#include <memory>
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

  string outputName;

  auto User_nEvts = INT_MIN;
  auto User_iEvt = INT_MIN;

  auto User_mem = INT_MIN;

  auto isVerbose = false;

  ProcessArgs(&theApp, &isVerbose,
			  &User_mem,
			  &User_nEvts, &User_iEvt,
			  &inputRATName, &inputFLATName,
			  &outputName);

  const string outName = outputName.empty() ?
						 ExtractFilenameFromPath(inputRATName) + "_WithRecon.root" : outputName;


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

  vector< pair<int, FlatRecon> > vpFlat
	  = GetVAssociatedIDAndRecon(inputFLATName.c_str(), isVerbose);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####      FILL RECON DATA TO ORIGINAL RAT MC           #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto tOutput = new TTree("OffT", "T With Offline Recon");
  auto BufDS = GetDS(inputRATName.c_str(), iEvt);
  AddReconInfo(BufDS,
			   0,
			   find_if(vpFlat.begin(), vpFlat.end(),
					   [&](const pair<int, FlatRecon>& p){return p.first == mMCID[iEvt++];})->second);
  tOutput->Branch("OffDS", "RAT::DS::Root", &BufDS);
  tOutput->Fill();

  OpenAndWriteTree(outName.c_str(), "RECREATE", tOutput);

  if(isVerbose)
	cout << "Writing recon info into DS object" << endl;

  ProgressBar progressBarRecon(nEvts, 70);

  unsigned long nBytes = 0;

  const unsigned long userSize = SetDefValue(User_mem, 100);
  const unsigned long MBSize = 1e6;
  const unsigned long bufferSize = userSize*MBSize;

  for(iEvt; iEvt<nEvts; iEvt++){

	// record the tick
	++progressBarRecon;

	int iBytes = 0;
	BufDS = GetDS(inputRATName.c_str(), iEvt, &iBytes);
	nBytes += iBytes;

	AddReconInfo(BufDS,
				 0,
				 find_if(vpFlat.begin(), vpFlat.end(),
						 [&](const pair<int, FlatRecon>& p){return p.first == mMCID[iEvt];})->second);

	tOutput->SetBranchAddress("OffDS", &BufDS, nullptr);
	tOutput->Fill();

	if(nBytes > bufferSize){

	  OpenAndWriteTree(outName.c_str(), "UPDATE", tOutput);

	  nBytes = 0;

	}

	if(isVerbose)
	  progressBarRecon.display();

  }

  if(isVerbose)
	progressBarRecon.done();

  OpenAndWriteTree(outName.c_str(), "UPDATE", tOutput);

  /////////////////////////
  // ...

  if(gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}

//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <csignal>
#include <numeric>
#include <climits>

/////////////////////////   ROOT   //////////////////////////
#include <TROOT.h>

/////////////////////////   USER   //////////////////////////
#include "utils.hh"
#include "FlattenHits.hh"

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
  string inputName;
  string outputName;

  auto User_nEvts = INT_MIN;
  auto User_iEvt = INT_MIN;

  auto isVerbose = false;
  auto isSource = false;

  ProcessArgs(&theApp, &isVerbose, &isSource,
			  &User_nEvts, &User_iEvt,
			  &inputName, &outputName);


  const string outName = outputName.empty() ?
						 ExtractFilenameFromPath(inputName) + "_FLAT.npz" : outputName;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE ANALYZER                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto FileAnalyzer = new Analyzer(inputName.c_str());


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                   ANALYSIS                        #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  unsigned long nEvts = User_nEvts > INT_MIN && User_nEvts < FileAnalyzer->GetNEvts() ?
						User_nEvts : FileAnalyzer->GetNEvts();
  unsigned long iEvt = SetDefValue(User_iEvt, 0);
  nEvts = iEvt > 0 ? iEvt + nEvts : nEvts;
  nEvts = nEvts > FileAnalyzer->GetNEvts() ? FileAnalyzer->GetNEvts() : nEvts;
  ProgressBar progressBar(nEvts, 70);

  GetVHitAndDumpFlatNPZ(FileAnalyzer, iEvt, outName, "w", isSource);
  iEvt++;

  for(iEvt; iEvt<nEvts; iEvt++) {

	if(EoF) break;

	// record the tick
	++progressBar;

	GetVHitAndDumpFlatNPZ(FileAnalyzer, iEvt, outName, "a", isSource);

	// display the bar
	if(isVerbose)
	  progressBar.display();

  }

  progressBar.done();

  cout << endl;
  EoF = 1;


  /////////////////////////
  // ...

  if(gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}

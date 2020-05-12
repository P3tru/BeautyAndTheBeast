//
// Created by zsoldos on 2/21/20.
//

#ifndef _FILLBACKRATMC_HH_
#define _FILLBACKRATMC_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>

/////////////////////////   BOOST  ///////////////////////////
#include <boost/algorithm/string/predicate.hpp>

/////////////////////////   ROOT   ///////////////////////////
#include <TApplication.h>
#include <TTree.h>

/////////////////////////   USER  ///////////////////////////
#include "utils.hh"

#include "Analyzer.hh"
#include "AnalyzerFunctions.hh"

#include "ProgressBar.hpp"

class OffEV : public TObject {

 public:
  OffEV() = default;
  virtual ~OffEV(){};

  RAT::DS::PathFit PF;
  int MCID;

};

typedef struct {double X, Y, Z, T, Theta, Phi;} FlatRecon;

TTree* SetFlatTreeReader(TFile *tFile,
						 FlatRecon &fRecon,
						 Long64_t &MCID){

  auto tTree = (TTree*)tFile->Get("Recon");

  tTree->SetBranchAddress("X", &fRecon.X);
  tTree->SetBranchAddress("Y", &fRecon.Y);
  tTree->SetBranchAddress("Z", &fRecon.Z);
  tTree->SetBranchAddress("T", &fRecon.T);
  tTree->SetBranchAddress("Theta", &fRecon.Theta);
  tTree->SetBranchAddress("Phi", &fRecon.Phi);
  tTree->SetBranchAddress("MCID", &MCID);

  return tTree;

};

void SetMC(RAT::DS::MC *MC, RAT::DS::MC *cMC){

  MC->SetID(cMC->GetID());
  MC->SetUTC(cMC->GetUTC());

  for(auto iPart=0; iPart<cMC->GetMCParticleCount(); iPart++){
    auto mcPart = MC->AddNewMCParticle();
    mcPart = dynamic_cast<RAT::DS::MCParticle *>(cMC->GetMCParticle(iPart)->Clone());
  }

  for(auto iPart=0; iPart<cMC->GetMCParentCount(); iPart++){
    auto mcPart = MC->AddNewMCParent();
    mcPart = dynamic_cast<RAT::DS::MCParticle *>(cMC->GetMCParent(iPart)->Clone());
  }

}

map<int, int> GetMAssociatedIDAndEvt(Analyzer *fA,
									 unsigned long iEvt, unsigned long nEvts, bool isVerbose){

  map<int, int> mMCID;

  if(isVerbose)
	cout << "Recovering DS address for each events" << endl;

  ProgressBar progressBar(nEvts, 70);

  for(iEvt; iEvt<nEvts; iEvt++) {

	// record the tick
	++progressBar;

	auto mc = GetRATMCOnEvt(fA, iEvt);
	mMCID[iEvt] = mc->GetID();

	if(isVerbose)
	  progressBar.display();

  }

  if(isVerbose)
	progressBar.done();

  return mMCID;

};

void AddReconInfo(RAT::DS::Root *DS, int iEV, FlatRecon fRecon){

  DS->GetEV(iEV)->GetPathFit()->SetTime(fRecon.T);
  DS->GetEV(iEV)->GetPathFit()->SetTime0(fRecon.T);
  DS->GetEV(iEV)->GetPathFit()->SetPos0(TVector3(fRecon.X, fRecon.Y, fRecon.Z));
  DS->GetEV(iEV)->GetPathFit()->SetPosition(TVector3(fRecon.X, fRecon.Y, fRecon.Z));
  TVector3 dir(1., 0. ,0.);
  dir.SetTheta(fRecon.Theta);
  dir.SetPhi(fRecon.Phi);
  dir.SetMag(1);
  DS->GetEV(iEV)->GetPathFit()->SetDirection(dir);

}

vector< pair<int, FlatRecon> > GetVAssociatedIDAndRecon(const char *fFLAT, bool isVerbose){

  vector< pair<int, FlatRecon> > vpFlat;

  auto tFile = TFile::Open(fFLAT);

  FlatRecon fRecon;
  Long64_t MCID;
  auto TRecon = SetFlatTreeReader(tFile,
								  fRecon,
								  MCID);

  auto nReconEntries = TRecon->GetEntries();

  for(auto iEntry=0; iEntry<nReconEntries; iEntry++) {

	TRecon->GetEntry(iEntry);
	vpFlat.emplace_back(make_pair(MCID, fRecon));

  }

  delete tFile;

  return vpFlat;

}

void ReadAndFillEvt(const char *fMC,
					unsigned long iEvt,
					unsigned nEV,
					FlatRecon vfRecon,
					RAT::DS::Root *DS){

  auto fA = new Analyzer(fMC);
  auto mc = GetRATMCOnEvt(fA, iEvt);
  DS = dynamic_cast<RAT::DS::Root *>(fA->GetDS()->Clone());

  delete fA;

  cout << DS->GetRunID() << endl;

  cout << vfRecon.X << " ";
  cout << vfRecon.Y << " ";
  cout << vfRecon.Z << endl;

  // for(auto iEV=0; iEV<nEV; iEV++){
  AddReconInfo(DS, 0, vfRecon);
  // }

  // tOutput->Fill();

}

RAT::DS::Root *GetDS(const char *fRATName, unsigned long iEvt, int *nBytes = nullptr){

  auto fMC = TFile::Open(fRATName, "READ");
  auto fMCTree = dynamic_cast<TTree *>(fMC->Get("T"));
  auto BufDS = new RAT::DS::Root();
  fMCTree->SetBranchAddress("ds", &BufDS);
  int nBytesBuf = fMCTree->GetEntry(iEvt);

  fMC->Close();
  delete fMC;

  if(nBytes)
    *nBytes = nBytesBuf;

  return BufDS;

}

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> -i FILE.root -o FILE.npz" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-v\tSet verbose mode\n"

	   << "\t-i--RAT\tinput RAT file (.root)\n"
	   << "\t-i--FLAT\tinput RECON file (.root)\n"

	   << "\t-NEvts\tNEvts to process (int)\n"
	   << "\t-iEvt\tStart at Evt #i (int)\n"

	   << endl;

}


void ProcessArgs(TApplication *theApp, bool *isVerbose,
				 int *User_nEvts, int *User_iEvt,
				 string *inputRATName, string *inputFLATName) {

  // Reading user input parameters
  if (theApp->Argc() < 2) {
	ShowUsage(theApp->Argv(0));
	exit(0);
  }

  for (int i = 1; i < theApp->Argc(); i++) {
	string arg = theApp->Argv(i);
	if ((arg == "-h") || (arg == "--help")) {
	  ShowUsage(theApp->Argv(0));
	  exit(0);

	} else if (boost::iequals(arg, "-v")) {
	  *isVerbose = true;

	} else if (boost::iequals(arg, "-NEvts")) {
	  *User_nEvts = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg, "-iEvt")) {
	  *User_iEvt = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg,"-i--RAT")) {
	  *inputRATName = theApp->Argv(++i);

	} else if (boost::iequals(arg,"-i--FLAT")) {
	  *inputFLATName = theApp->Argv(++i);

	} else {
	  cout << "Unkown parameter" << endl;
	  continue;
	}
  }

  if(inputRATName->empty() && inputFLATName->empty()){
	cout << "ERROR: No input file provided!" << endl;
	exit(EXIT_FAILURE);
  } else if(!IsFileExist(*inputRATName) && !IsFileExist(*inputFLATName)){
	cout << "ERROR: file doesn't exist!!" << endl;
	exit(EXIT_FAILURE);
  }

}

#endif //_FILLBACKRATMC_HH_

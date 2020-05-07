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

class OffEV : public TObject {

 public:
  OffEV() = default;
  virtual ~OffEV(){};

  RAT::DS::PathFit PF;
  int MCID;

};


TTree* SetFlatTreeReader(TFile *tFile,
						 double &X, double &Y, double &Z,
						 double &T,
						 double &Theta, double &Phi,
						 Long64_t &MCID){

  auto tTree = (TTree*)tFile->Get("Recon");

  tTree->SetBranchAddress("X", &X);
  tTree->SetBranchAddress("Y", &Y);
  tTree->SetBranchAddress("Z", &Z);
  tTree->SetBranchAddress("T", &T);
  tTree->SetBranchAddress("Theta", &Theta);
  tTree->SetBranchAddress("Phi", &Phi);
  tTree->SetBranchAddress("MCID", &MCID);

  return tTree;

};

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

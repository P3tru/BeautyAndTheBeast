//
// Created by zsoldos on 2/21/20.
//

#ifndef _FLATTENHITS_HH_
#define _FLATTENHITS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>

/////////////////////////   BOOST  ///////////////////////////
#include <boost/algorithm/string/predicate.hpp>

/////////////////////////   ROOT   ///////////////////////////
#include <TApplication.h>

/////////////////////////   USER  ///////////////////////////
#include "utils.hh"

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> -i FILE.root -o FILE.npz" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-v\tSet verbose mode\n"

	   << "\t-t\tSet if tracks were recorded\n"

	   << "\t-i\tinput  file (.root)\n"
	   << "\t-o\toutput file (.npz)\n"

	   << "\t-NEvts\tNEvts to process (int)\n"
	   << "\t-iEvt\tStart at Evt #i (int)\n"

	   << endl;

}


void ProcessArgs(TApplication *theApp, bool *isVerbose, bool *isSource,
				 int *User_nEvts, int *User_iEvt,
				 string *filename, string *outputname) {

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

	} else if (boost::iequals(arg, "-t")) {
	  *isSource = true;

	} else if (boost::iequals(arg, "-NEvts")) {
	  *User_nEvts = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg, "-iEvt")) {
	  *User_iEvt = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg,"-i")) {
	  *filename = theApp->Argv(++i);

	} else if (boost::iequals(arg,"-o")) {
	  *outputname = theApp->Argv(++i);

	} else {
	  cout << "Unkown parameter" << endl;
	  continue;
	}
  }

  if(filename->empty()){
	cout << "ERROR: No input file provided!" << endl;
	exit(EXIT_FAILURE);
  } else if(!IsFileExist(*filename)){
	cout << "ERROR: file doesn't exist!!" << endl;
	exit(EXIT_FAILURE);
  }

}

#endif //_FLATTENHITS_HH_

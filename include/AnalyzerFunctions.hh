//
// Created by zsoldos on 2/21/20.
//

#ifndef _ANALYZERFUNCTIONS_HH_
#define _ANALYZERFUNCTIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   USER   //////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"
#include "FlatParticle.hh"

using namespace std;

// Add list of analyzer to vector
void AddFAnalyzers(vector<Analyzer*> *vFAnalyzer, const string& inputName, const string& listName = "");

// Access the MC object from RAT
// at iEvt (from MCTree)
RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0);

// Fill vector<Hit> with MC Info
vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt=0, bool isSource=false);

// Fill vector<Hit> with MC Info from one mother particle
vector<Hit> GetVHitsFromPart(Analyzer *fAnalyzer, unsigned int iEvt=0,
							 string sPartName="e+",
							 vector<ComplexParticle> *vCPart = NULL);

// Dump vector<Hit> inside a npz file
void GetVHitAndDumpFlatNPZ(Analyzer *fAnalyzer, unsigned iEvt, const string& NPZName, const string& mode="a");

#endif
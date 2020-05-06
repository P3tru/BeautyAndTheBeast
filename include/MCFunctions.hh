//
// Created by zsoldos on 2/24/20.
//

#ifndef _MCFUNCTIONS_HH_
#define _MCFUNCTIONS_HH_


///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>

/////////////////////////   ROOT  ///////////////////////////
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"
#include "FlatParticle.hh"

using namespace std;

// Print basic track info
// Useful for debugging
void PrintTrackInfo(RAT::DS::MC *mc, unsigned int iTrack);

// Get NumPE and NHits from MC object
void GetMCNHitsAndNPE(double *NHits, double *NPE, RAT::DS::MC * mc);

// Get vector<FlatPhoton> of photons in evt, storing their wl and create proc
vector<FlatPhoton> GetPhotonsFromEvt(RAT::DS::MC * mc);

// Get info from MC particle in a flat way
void GetPartInfoFromTrackStep(FlatParticle *fp, RAT::DS::MCTrack *mctrack);

// Get ID from MC particle in a flat way
void GetPartIDFromTrackStep(G4Particle *p, RAT::DS::MCTrack *mctrack);

// Get vector<ComplexParticle> of MC object
vector<ComplexParticle> GetVPart(RAT::DS::MC * mc, const string& sPartName);

// Fill vector<ComplexParticle> with daughters photons
vector<FlatPhoton> GetAndFillPhotonsToVPart(RAT::DS::MC * mc, vector<ComplexParticle> *vPart);

bool IsParticle(RAT::DS::MCTrack *mctrack, const string& name);

#endif
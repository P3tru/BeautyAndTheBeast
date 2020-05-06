//
// Created by zsoldos on 2/24/20.
//

#ifndef _HITFUNCTIONS_HH_
#define _HITFUNCTIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   ROOT  ///////////////////////////
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "HitClass.hh"

using namespace std;

// Operators to add and substract hit
// ONLY used for shifting the hits in TIME.
// which is useful for events where RAT doesn't split them
Hit operator-(Hit h1, Hit h2);
Hit operator+(Hit h1, Hit h2);

// Print all hits information: Pos, Q, T
void PrintVHits(vector<Hit> const Hits);

// Sort hits from Time
void SortVHits(vector<Hit> *vHit);

// Split long hits collection
// If a collection of hits is superior than DAQWindow,
// return a new vector<Hit> and erase from vHit the hits after DAQWindow
vector<Hit> SplitVHits(vector<Hit> *vHit, double DAQWindow);

// Return a new vector<Hit> where time is corrected from the first hit,
// and use can add a pretrig time (default time=0.)
vector<Hit> ResetTVHits(vector<Hit> rawHits,
						const Hit& PreTrig = Hit(TVector3(0.,0.,0.), 0., 0.));

// Erase from vector<Hit> hits after hCut
void RemoveHitsAfterCut(vector<Hit> &Hits,
						const Hit& hCut = Hit(TVector3(0.,0.,0.), 0., 0.));

// Fill residuals from vector<Hit>
// 1D or 2D, with the CosTheta of the hits with respect to particle direction
void FillResiduals(TH1D *hResid, const vector<Hit>& vHit,
				   const TVector3& Origin = TVector3(0.,0.,0.),
				   double SoL = 224.9,
				   bool isWeight=false);
void FillCTheta(TH1D *hCT, const vector<Hit>& vHit,
				const TVector3& TrueDir = TVector3(0.,0.,0.));
void FillTResVSCosTheta(TH2D *h2D, const vector<Hit>& vHit,
						const TVector3& TrueOrigin = TVector3(0.,0.,0.),
						double SoL = 224.9,
						const TVector3& TrueDir = TVector3(0.,0.,0.),
						bool isWeight = false);
void FillAll(TH1D *hResid, TH1D *hCT, TH2D *h2D, const vector<Hit>& vHit,
			 const TVector3& TrueOrigin = TVector3(0.,0.,0.),
			 double SoL = 224.9,
			 const TVector3& TrueDir = TVector3(0.,0.,0.),
			 bool isWeight=false);

// Create histogram Q VS T
TH2D *GetHQVST(vector<Hit> vHit);

// Get directly from vector<Hit> NPE and NHits
void GetNPEAndNHitsFromHits(vector<Hit> Hits, double *NPE, double *NHits);

// Flatten hit vector for numpy use
vector<double> FlatenVHit(const vector<Hit>& vHit, int MCID=-1);

#endif //_HITFUNCTIONS_HH_
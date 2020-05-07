///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>
#include <fstream>

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "AnalyzerFunctions.hh"
#include "HitClass.hh"
#include "HitFunctions.hh"
#include "MCFunctions.hh"

#include "cnpy.h"

void AddFAnalyzers(vector<Analyzer*> *vFAnalyzer, const string& inputName, const string& listName){

  if(!inputName.empty()){
	vFAnalyzer->push_back(new Analyzer(inputName.c_str()));
  }

  if(!listName.empty()) {

	// Get ready to read inside inputName
	// each file name
	string line;
	ifstream file(listName);

	while (getline(file, line)) {

	  if (line.compare(0, 1, "#") == 0) {

		continue;

	  } else {

		vFAnalyzer->push_back(new Analyzer(line.c_str()));

	  }

	}

  }


}

RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetMC();

}

RAT::DS::EV * GetRATEVOnEvt(Analyzer *fAnalyzer, unsigned int iEvt, unsigned iEV){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetEV(iEV);

}

vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt, bool isSource){

  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  vector<Hit> vHit;

  vector<FlatPhoton> vFP;
  if(isSource){
	vFP = GetPhotonsFromEvt(mc);
  }

  for (int iPMT = 0; iPMT < mc->GetMCPMTCount(); iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  auto isDarkHit = mcPhoton->IsDarkHit();
	  if(isDarkHit)
		cout << "DARK HIT" << endl;

	  if(isSource){

		auto trackID = mcPhoton->GetTrackID();

		for(auto p:vFP){

		  if(trackID == p.GetTrackId()){

			// ########## //
			// FILL EVENT //
			// ########## //

			auto CreaProc = p.GetProcess();
			auto WL = p.GetWl();

			vHit.emplace_back(Hit(PMTPos, Q, T,
								  TVector3(0.,0.,0.),
								  TVector3(0.,0.,0.),
								  CreaProc, WL));

			break;

		  }

		}

	  } else {

		// ########## //
		// FILL EVENT //
		// ########## //

		vHit.emplace_back(Hit(PMTPos, Q, T));

	  }

	}

  } // END FOR iPMT

  return vHit;

}

vector<Hit> GetVHitsFromPart(Analyzer *fAnalyzer, unsigned int iEvt,
							 string sPartName,
							 vector<ComplexParticle> *vCPart){

  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  vector<ComplexParticle> vPart = GetVPart(mc, sPartName);
  vector<FlatPhoton> vFP = GetAndFillPhotonsToVPart(mc, &vPart);

  vector<Hit> vHit;

  for (int iPMT = 0; iPMT < mc->GetMCPMTCount(); iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  // cout << mcPhoton->GetTrackID() << endl;

	  // ########## //
	  // FILL EVENT //
	  // ########## //

	  bool isPhotonProcessed = false;

	  for (auto proton = vPart.begin(); proton != vPart.end(); proton++) {

		for (auto photon = proton->GetVp().begin(); photon != proton->GetVp().end(); photon++) {

		  if (mcPhoton->GetTrackID() == photon->GetTrackId()) {

			vHit.emplace_back(Hit(PMTPos, Q, T, proton->GetPos(), proton->GetDir()));

			isPhotonProcessed = true;
			proton->RemovePhoton(photon);
			break;

		  }

		  if (isPhotonProcessed)
			break;

		} // END for photon

	  } // END for proton

	} // END FOR iPhoton

  } // END for iPMT

  if(vCPart)
	*vCPart = vPart;

  return vHit;

}

void GetVHitAndDumpFlatNPZ(Analyzer *fAnalyzer, unsigned iEvt, const string& NPZName, const string& mode){

  vector<Hit> vHit = GetMCHitCollection(fAnalyzer, iEvt, true);
  auto MCID = GetRATMCOnEvt(fAnalyzer, iEvt)->GetID();

  const auto NHits = vHit.size();
  vector<double> vNPY = FlatenVHit(vHit, MCID);

  cnpy::npz_save(NPZName,
				 Form("Evt%d", iEvt),
				 &vNPY[0],
				 {NHits, 6}, // NHits vector of 6 dimension {X, Y, Z, T, Source, MCID}
				 mode);

}

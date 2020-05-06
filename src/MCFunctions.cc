//
// Created by zsoldos on 2/21/20.
//

///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   ROOT  ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "MCFunctions.hh"

using namespace std;

// Convert photon MeV to nm
static double MeV2lambda(double MeV){

  const double hc = 1.23984193e-3; //MeV.nm
  return hc/MeV;

}

bool operator==(G4Particle const& a, G4Particle const& b){
  return a.GetTrackId() == b.GetTrackId() && a.GetParentId() == b.GetParentId();
}

static int ConvertProcess(string proc){
  if(proc == "Cerenkov"){
	return 0;
  } else if (proc == "Scintillation") {
	return 1;
  } else if(proc.rfind("Reemission", 0) == 0) {
	return 2;
  } else {
	return -1;
  }
}

vector<FlatPhoton> GetPhotonsFromEvt(RAT::DS::MC * mc){

  vector<FlatPhoton> vFP;

  for(auto iTrack=0; iTrack<mc->GetMCTrackCount(); iTrack++){

    auto track = mc->GetMCTrack(iTrack);

    if(IsParticle(track, "opticalphoton")){

	  auto first = track->GetMCTrackStep(0);
	  FlatPhoton p;
	  GetPartIDFromTrackStep(&p, track);
	  p.SetWl(MeV2lambda(first->GetKE()));
	  p.SetProcess(ConvertProcess(first->GetProcess()));

	  vFP.emplace_back(p);

    }

  }

  return vFP;

}

void GetPartInfoFromTrackStep(FlatParticle *fp, RAT::DS::MCTrack *mctrack){

  auto *mctrackstep = mctrack->GetMCTrackStep(0);

  fp->SetPos(mctrackstep->GetEndpoint());
  fp->SetKinE(mctrackstep->GetKE());
  fp->SetDir(mctrackstep->GetMomentum().Unit());

}

void GetPartIDFromTrackStep(G4Particle *p, RAT::DS::MCTrack *mctrack){

  p->SetParentId(mctrack->GetParentID());
  p->SetTrackId(mctrack->GetID());

}

vector<ComplexParticle> GetVPart(RAT::DS::MC * mc, const string& sPartName){

  auto NbTracks = mc->GetMCTrackCount();

  vector<ComplexParticle> vPart;

  // First, identify how many particle with name sPartName are created in the evt.
  // Loop over tracks and records the one called sPartName.
  // Save only information of interest.

  for(auto iTrack=0; iTrack<NbTracks; iTrack++) {

	auto *mctrack = mc->GetMCTrack(iTrack);

	if (IsParticle(mctrack, sPartName)) {

	  ComplexParticle cPart;
	  cPart.SetName(mctrack->GetName());
	  GetPartIDFromTrackStep(&cPart, mctrack);
	  GetPartInfoFromTrackStep(&cPart, mctrack);

	  vPart.emplace_back(cPart);

	}

  }

  return vPart;

}

vector<FlatPhoton> GetAndFillPhotonsToVPart(RAT::DS::MC * mc, vector<ComplexParticle> *vPart){

  vector<FlatPhoton> vFP;

  auto NbTracks = mc->GetMCTrackCount();

  for(auto iTrack=0; iTrack<NbTracks; iTrack++) {

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "opticalphoton")) {

	  FlatPhoton fp;
	  GetPartIDFromTrackStep(&fp, mctrack);
	  fp.SetWl(MeV2lambda(mctrack->GetMCTrackStep(0)->GetKE()));

	  for(auto &p:*vPart){

		if(fp.GetParentId() == p.GetTrackId()){

		  p.AddPhoton(fp);
		  vFP.emplace_back(fp);

		} // END if daughter

	  } // ENF for vPart

	} // END if photon

  } // END for iTrack

  return vFP;

}


bool IsParticle(RAT::DS::MCTrack *mctrack, const string& name){

  return mctrack->GetParticleName() == name;

}


void PrintTrackInfo(RAT::DS::MC *mc, unsigned int iTrack){

  auto * mctrack = mc->GetMCTrack(iTrack);
  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

  cout << mctrack->GetParticleName() << endl;
  cout << " With ID: " << mctrack->GetID() << " and parent ID: " << mctrack->GetParentID() << endl;
  cout << "TOT Nb Steps: " << NbTrackSteps << endl;

  for(auto iStep=0; iStep<NbTrackSteps; iStep++){

	auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

	cout << mctrackstep->GetLength() <<"mm " << mctrackstep->GetKE() << "MeV " << endl;

  }


}

void GetMCNHitsAndNPE(double *NHits, double *NPE, RAT::DS::MC * mc){

  *NHits=mc->GetMCPMTCount();
  *NPE=mc->GetNumPE();

}


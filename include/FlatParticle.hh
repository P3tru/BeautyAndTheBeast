//
// Created by zsoldos on 4/10/20.
//

#ifndef _FLATPARTICLE_HH_
#define _FLATPARTICLE_HH_

class G4Particle{
 protected:
  int parentID;
  int trackID;
 public:
  G4Particle() = default;
  G4Particle(int parent_id, int track_id) : parentID(parent_id), trackID(track_id) {}
  virtual ~G4Particle() {}
  int GetParentId() const { return parentID; }
  void SetParentId(int parent_id) { parentID = parent_id; }
  int GetTrackId() const { return trackID; }
  void SetTrackId(int track_id) { trackID = track_id; }
};

class FlatParticle: public G4Particle{
 protected:
  string name;
  TVector3 pos;
  TVector3 dir;
  double KinE;
 public:
  FlatParticle() = default;
  FlatParticle(const string &name, const TVector3 &pos, const TVector3 &dir, double kin_e)
	  : name(name), pos(pos), dir(dir), KinE(kin_e) {
	parentID=-1;
	trackID=-1;
  }
  FlatParticle(int parent_id, int track_id, const string &name, const TVector3 &pos, const TVector3 &dir, double kin_e)
	  : G4Particle(parent_id, track_id), name(name), pos(pos), dir(dir), KinE(kin_e) {}
  virtual ~FlatParticle() {}
  const string &GetName() const { return name; }
  void SetName(const string &name) { FlatParticle::name = name; }
  const TVector3 &GetPos() const { return pos; }
  void SetPos(const TVector3 &pos) { FlatParticle::pos = pos; }
  const TVector3 &GetDir() const { return dir; }
  void SetDir(const TVector3 &dir) { FlatParticle::dir = dir; }
  double GetKinE() const { return KinE; }
  void SetKinE(double kin_e) { KinE = kin_e; }
};

class FlatPhoton: public G4Particle{
 protected:
  double wl;
  int process;
 public:
  FlatPhoton() = default;
  FlatPhoton(int parent_id, int track_id, double wl)
	  : G4Particle(parent_id, track_id), wl(wl), process(-1) {}
  FlatPhoton(int parent_id, int track_id, double wl, int process)
	  : G4Particle(parent_id, track_id),
		wl(wl),
		process(process) {}
  virtual ~FlatPhoton() {}
  double GetWl() const { return wl; }
  void SetWl(double wl) { FlatPhoton::wl = wl; }
  int GetProcess() const { return process; }
  void SetProcess(int process) { FlatPhoton::process = process; }
};

class ComplexParticle: public FlatParticle{
 protected:
  double dE;
  double dX;
  vector<FlatPhoton> vP;
 public:
  ComplexParticle() = default;
  ComplexParticle(const string &name, const TVector3 &pos, const TVector3 &dir, double kin_e)
	  : FlatParticle(name, pos, dir, kin_e) {
    dE=0;
	dX=0.;
  }
  ComplexParticle(int parent_id, int track_id,
				  const string &name,
				  const TVector3 &pos, const TVector3 &dir,
				  double kin_e)
	  : FlatParticle(parent_id, track_id, name, pos, dir, kin_e) {
	dE=0;
	dX=0.;
  }
  virtual ~ComplexParticle() { }
  void AddPhoton(FlatPhoton p){
	vP.emplace_back(p);
  }
  void RemovePhoton(__gnu_cxx::__normal_iterator<const FlatPhoton *, vector<FlatPhoton>> itPhoton){
    vP.erase(itPhoton);
  }

  const vector<FlatPhoton> &GetVp() const { return vP; }
  double GetDe() const { return dE;}
  void SetDe(double d_e) { dE = d_e; }
  double GetDx() const { return dX; }
  void SetDx(double d_x) { dX = d_x; }
  void ELoss(double E){ dE+=E; }
  void AddTrackLength(double X) { dX=+X; }

};


#endif //_FLATPARTICLE_HH_

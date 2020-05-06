//
// Created by zsoldos on 2/24/20.
//

#ifndef _HITCLASS_HH_
#define _HITCLASS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include <TVector3.h>


#define PI 3.14159
#define C 299.792458 // mm/ns

using namespace std;

class THit {

};

class Hit {

 protected:
  TVector3 Pos;

  double Q;
  double T;

  double TResid;
  double D;

  TVector3 TrueOrigin;
  TVector3 TrueDir;

  int TrueProcess;
  double TrueWL;

 public:
  Hit()
	  : Pos(TVector3(0.,0.,0.)) ,
		Q(-1), T(-1),
		TResid(-1), D(-1),
		TrueOrigin(TVector3(0.,0.,0.)), TrueDir(TVector3(0.,0.,0.)),
		TrueProcess(-1), TrueWL(-1.) { }
  Hit(const TVector3 &pos)
	  : Pos(pos),
		Q(-1), T(-1),
		TResid(-1), D(-1),
		TrueOrigin(TVector3(0.,0.,0.)), TrueDir(TVector3(0.,0.,0.)),
		TrueProcess(-1), TrueWL(-1.) { }
  Hit(const TVector3 &pos, double q, double t)
	  : Pos(pos),
		Q(q), T(t),
		TResid(-1), D(-1),
		TrueOrigin(TVector3(0.,0.,0.)), TrueDir(TVector3(0.,0.,0.)),
		TrueProcess(-1), TrueWL(-1.) { }
  Hit(const TVector3 &pos, double q, double t, int true_process, double true_wl)
	  : Pos(pos),
		Q(q), T(t),
		TResid(-1), D(-1),
		TrueOrigin(TVector3(0.,0.,0.)), TrueDir(TVector3(0.,0.,0.)),
		TrueProcess(true_process), TrueWL(true_wl) { }
  Hit(const TVector3 &pos, double q, double t, const TVector3 &true_origin, const TVector3 &true_dir)
	  : Pos(pos), Q(q), T(t),
		TResid(-1), D(-1),
		TrueOrigin(true_origin), TrueDir(true_dir),
		TrueProcess(-1), TrueWL(-1.) { }
  Hit(const TVector3 &pos,
	  double q, double t,
	  const TVector3 &true_origin, const TVector3 &true_dir,
	  int true_process, double true_wl)
	  : Pos(pos),
	  Q(q), T(t),
	  TrueOrigin(true_origin), TrueDir(true_dir),
	  TrueProcess(true_process), TrueWL(true_wl) {}

  const TVector3 &GetPos() const { return Pos; }
  void SetPos(const TVector3 &pos) { Pos = pos; }

  double GetQ(){ return Q; };
  void SetQ(double x){ Q = x; };
  double GetT(){ return T; };
  void SetT(double x){ T = x; };

  double GetTResid() const { return TResid; }
  void SetTResid(double t_resid) { TResid = t_resid; }
  double GetD() const { return D; }
  void SetD(double d) { D = d; }

  const TVector3 &GetTrueOrigin() const { return TrueOrigin; }
  void SetTrueOrigin(const TVector3 &true_origin) { TrueOrigin = true_origin; }
  const TVector3 &GetTrueDir() const { return TrueDir; }
  void SetTrueDir(const TVector3 &true_dir) { TrueDir = true_dir; }

  int GetTrueProcess() const { return TrueProcess; }
  void SetTrueProcess(int true_process) { TrueProcess = true_process; }
  double GetTrueWl() const { return TrueWL; }
  void SetTrueWl(double true_wl) { TrueWL = true_wl; }

  double CalculateDistance(const TVector3& Origin = TVector3(0,0,0)){
	return TVector3(Pos-Origin).Mag();
  };
  double CalculateTResid(const TVector3& Origin = TVector3(0,0,0), double SoL = 299.792458){
	return T - CalculateDistance(Origin)/SoL;
  };

  double CalcuateCosTHit(const TVector3& Origin = TVector3(0,0,0),
						 const TVector3& Dir = TVector3(0,0,0)){
    return cos(Dir.Angle(Pos-Origin));
  }

  bool operator==(const Hit &rhs) const { return T == rhs.T; }
  bool operator!=(const Hit &rhs) const { return !(rhs == *this); }
  bool operator<(const Hit &rhs) const { return T < rhs.T; }
  bool operator>(const Hit &rhs) const { return rhs < *this; }
  bool operator<=(const Hit &rhs) const { return !(rhs < *this); }
  bool operator>=(const Hit &rhs) const { return !(*this < rhs); }

  void operator+=(Hit rhs) { SetT(T+rhs.GetT()); }
  void operator-=(Hit rhs) { SetT(T-rhs.GetT()); }

  void Print() const {
	cout << "X: " << Pos.x()
		 << " Y: " << Pos.y()
		 << " Z: " << Pos.z()
		 << " Q: " << Q
		 << " T: " << T << endl;
  }
  void PrintTrue() const {
    TrueOrigin.Print();
    TrueDir.Print();
  }

};

// Applied only on time; Super useful to do cuts.
Hit operator-(Hit h1, Hit h2);
Hit operator+(Hit h1, Hit h2);

class HitCut{

 public:
  Hit hCut;

  explicit HitCut(const Hit &h_cut) : hCut(h_cut) {}

  int operator()(Hit h){
	bool YesItIs = false;
	if(h>hCut){
	  YesItIs = true;
	}
	return YesItIs;
  }

};


#endif //_HITCLASS_HH_

//
// Created by zsoldos on 2/21/20.
//

#ifndef _ANALYZER_HH_
#define _ANALYZER_HH_

/////////////////////////   ROOT   //////////////////////////
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>

using namespace std;

class Analyzer {

 protected:
  TFile *fMC;
  TTree *treeMC;
  RAT::DS::Root *rds;

  RAT::DS::Run *run;
  TTree *runTree;

  unsigned long nEvts;

  string tag;

 public:
  Analyzer() = default;

  Analyzer(TFile *f_mc,
		   TTree *tree_mc,
		   RAT::DS::Root *rds,
		   RAT::DS::Run *run,
		   TTree *run_tree,
		   unsigned long n_evts,
		   const string &tag)
	  : fMC(f_mc), treeMC(tree_mc), rds(rds), run(run), runTree(run_tree), nEvts(n_evts), tag(tag) {

  }

  explicit Analyzer(const char* filename, const string &tag=""){

	if(tag.empty())
	  Analyzer::tag.clear();
	else
	  Analyzer::tag = tag;

	SetAnalyzer(filename);

  }

  void SetAnalyzer(const char *filename){

	Analyzer::fMC = new TFile(filename, "UPDATE");

	if (Analyzer::fMC->IsOpen()) {

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  Analyzer::run = new RAT::DS::Run();

	  Analyzer::runTree = (TTree *) Analyzer::fMC->Get("runT");
	  Analyzer::runTree->SetBranchAddress("run", &run);
	  Analyzer::runTree->GetEntry(0);

	  Analyzer::treeMC = (TTree *) Analyzer::fMC->Get("T");

	  Analyzer::rds = new RAT::DS::Root();
	  Analyzer::treeMC->SetBranchAddress("ds", &rds);

	  Analyzer::nEvts = Analyzer::treeMC->GetEntries();

	}

  }

  virtual ~Analyzer() {
	delete rds;
	delete treeMC;

	delete runTree;
	delete run;

	fMC->Close();
	delete fMC;

  }

  TFile *GetFmc() const {
	return fMC;
  }
  TTree *GetTreeMc() const {
	return treeMC;
  }
  RAT::DS::Root *GetRds() const {
	return rds;
  }

  unsigned long GetNEvts() const {
	return nEvts;
  }

  RAT::DS::Run *GetRun() const {
	return run;
  }
  TTree *GetRunTree() const {
	return runTree;
  }

  const string &GetTag() const {
	return tag;
  }
  void SetTag(const string &tag) {
	Analyzer::tag = tag;
  }

};

#endif //
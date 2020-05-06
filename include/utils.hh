//
// Created by zsoldos on 11/26/19.
//

#ifndef _UTILS_HH_
#define _UTILS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TStyle.h>
#include <TVector3.h>
#include <TH2D.h>
#include <TFile.h>

/////////////////////////   RAT   ///////////////////////////

using namespace std;

#define SQRT2 1.41421356237

static string ExtractFilenameFromPath(string pathname){

  boost::filesystem::path p(pathname);
  return p.filename().string();

}

int EoF = 0;

static void Interrupt(int arg){

  if(EoF==0) { printf("got a control-C, stop\n"); EoF=1; return; }
  else { printf("got a control-C, exiting\n"); exit(0); }

}

static void SetBasicStyle(){
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(0); // This determines if you want a stats box
  gStyle->SetOptFit(0); // This determines if you want a fit info box
  gStyle->GetAttDate()->SetTextColor(1);
  gStyle->SetOptTitle(1); // no title; comment out if you want a title
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132,"XYZ");

  gROOT->ForceStyle();

  gStyle->SetPalette(kDarkRainBow);
}

template <typename T>
static void SetBasicTStyle(T *h,
			   int Color=kBlue-4,
			   int LineWidth=1, int LineStyle=1,
			   int MarkerSize=1, int MarkerStyle=kPlus){

  h->SetLineColor(Color);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);

  h->SetMarkerColor(Color);
  h->SetMarkerSize(MarkerSize);
  h->SetMarkerStyle(MarkerStyle);

}

static bool IsFileExist(const string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

template <typename T>
T SetDefValue(const T& User, const T& def){
  return (User > std::numeric_limits<T>::min() ) ? User : def;
}

// define a function with 3 parameters
static Double_t fitGaus(Double_t *x,Double_t *par) {

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

// define a function with 3 parameters
static Double_t fitExpo(Double_t *x,Double_t *par) {

  return par[0]*TMath::Exp(-x[0]/par[1]);

}

static Double_t fitGausExpoTail(Double_t *x,Double_t *par){

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  if(x[0]>0)
    fitval += par[3]*TMath::Exp(-x[0]/par[4]);
  return fitval;

}

static Double_t fitEMG(Double_t *x,Double_t *par){

  Double_t xx = x[0];
  Double_t norm = par[0];
  Double_t mu = par[1];
  Double_t sigma = par[2];
  Double_t lambda = par[3];

  Double_t argNorm = norm*lambda/2;
  
  Double_t argExpo = (lambda/2)*(2*mu + lambda*sigma*sigma -2*xx);

  Double_t argErf = (mu + lambda*sigma*sigma - xx)/(SQRT2*sigma);

  return argNorm * TMath::Exp(argExpo) * TMath::Erfc(argErf);

}

template <typename T>
TCanvas *PlotAHist(T *h, const char *opt=""){

  auto *c1 = new TCanvas(Form("c%s", h->GetName()), Form("c%s", h->GetName()), 800, 600);
  c1->SetGrid();
  h->Draw(opt);
  return c1;

}

template <typename T>
double CalculateProb(T *hPDF, T *hExp, double *Chi2 = NULL, int *NdF = NULL){

  if(!hPDF || !hExp)
    return -1;

  if(hPDF->GetEntries() == 0 || hExp->GetEntries() == 0)
    return -1;

  auto nBinsX = hPDF->GetNbinsX();
  auto nBinsY = hPDF->GetNbinsY();

  auto N = hExp->Integral();
  auto W = hPDF->GetSumOfWeights();
  auto normW = 2*W*W;

  double chi2 = 0.;
  int NonNullBin = 0;

  for(auto iBinX=1; iBinX<=nBinsX; iBinX++){
    for(auto iBinY=1; iBinY<=nBinsY; iBinY++) {

      double n = hExp->GetBinContent(hExp->GetBin(iBinX, iBinY));
      double w = hPDF->GetBinContent(hPDF->GetBin(iBinX, iBinY));

      if(n == 0 || w == 0) continue;

      double s2 = pow(hPDF->GetBinError(hPDF->GetBin(iBinX, iBinY)),2);
      double res = W*w - N*s2;

      double P = res + sqrt(pow(res, 2) + 4*W*W*s2*s2*n);
      if(P == 0) continue;
      P/=normW;

      chi2+=pow(n - N*P,2)/(N*P) + pow(w - W*P,2)/s2;
      NonNullBin++;

    }
  }

  if(Chi2)
    *Chi2=chi2;
  if(NdF)
    *NdF=NonNullBin-1;

  return TMath::Prob(chi2, NonNullBin-1);


}

#endif // _UTILS_HH_

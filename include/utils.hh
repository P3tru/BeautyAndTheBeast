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
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>

/////////////////////////   RAT   ///////////////////////////

using namespace std;

#define SQRT2 1.41421356237

// void loadlibs(){

//   R__LOAD_LIBRARY(libEVFunctions.so);
//   R__LOAD_LIBRARY(libMCFunctions.so);
//   R__LOAD_LIBRARY(libHitFunctions.so);
//   R__LOAD_LIBRARY(libLL.so);
//   R__LOAD_LIBRARY(libCalibFunctions.so);
//   R__LOAD_LIBRARY(libAnalyzerFunctions.so);
//   R__LOAD_LIBRARY(libcnpy.so);

// }

static string ExtractFilenameFromPath(string pathname){

  boost::filesystem::path p(pathname);
  return p.filename().string();

}

static vector<double> CorrectBinRangeArray(vector<double> inputArray){

  vector<double> output;
  output.resize(inputArray.size()+1);

  double DeltaBin0 = inputArray[1] - inputArray[0];

  output[0] = inputArray[0] - DeltaBin0/2;

  for(int i = 1; i<inputArray.size(); i++){

    output[i] = inputArray[i-1] + (inputArray[i] - inputArray[i-1])/2;

  }

  output[inputArray.size()] =
    inputArray[inputArray.size()-1] +
    (inputArray[inputArray.size()-1] - inputArray[inputArray.size()-2])/2;

  return output;

}

int EoF = 0;

static void Interrupt(int arg){

  if(EoF==0) { printf("got a control-C, stop\n"); EoF=1; return; }
  else { printf("got a control-C, exiting\n"); exit(0); }

}

TH1D *GetHPDF(const char *filename, const char *histname) {

  TFile *fPDF = new TFile(filename);
  TH1D *hPDF = nullptr;

  if(fPDF){

    hPDF = (TH1D *) fPDF->Get(histname)->Clone();

  }

  return hPDF;

}

static void SetStyleVariables(TStyle *t2kStyle){

  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.15);
  t2kStyle->SetPadRightMargin(0.15); //0.05
  t2kStyle->SetPadBottomMargin(0.16);
  t2kStyle->SetPadLeftMargin(0.13);

  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.05);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  // t2kStyle->SetPadGridX(true);
  // t2kStyle->SetPadGridY(true);


  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes


  // t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs,
				   stops, red, green, blue,
				   NCont);
  t2kStyle->SetNumberContours(NCont);

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

static void SetBasicTH1Style(TH1 *h,
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

double GetBetaFromEKin(double EKin, double Mass){
  return TMath::Sqrt(EKin * (EKin + 2*Mass)) / (EKin + Mass);
}

double GetBetheBloch(double Beta){

  const double re = 2.817940325; //fm
  const double mass_e = 0.510998918; //c2 => Mev
  const double K = 4 * TMath::Pi() * TMath::Na() * re*re * mass_e*mass_e;
  const double Z = 0;
  const double A = 0;
  const double I = 0;
  const double TMax = 0;
  const double z = 1;

  const double gamma = 1 / TMath::Sqrt(1 - Beta*Beta);

  return K * z*z * (Z/A) * (1/(Beta*Beta)) * (0.5*TMath::Log(2*mass_e*Beta*Beta*gamma*gamma*TMax/(I*I)) - Beta*Beta);

}

#endif // _UTILS_HH_

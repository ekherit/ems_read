/*
 * =====================================================================================
 *
 *       Filename:  read_ems.cpp
 *
 *    Description:  Read E.results  and P.results
 *
 *        Version:  1.0
 *        Created:  25.04.2018 16:22:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <list>
#include <fmt/printf.h>
#include <string>
#include <limits>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TF1.h>
#include <TSpline.h>
#include <TSystem.h>


#include <ibn/valer.h>

#include "multistep.h"

struct EMS_t
{
  ibn::valer<time_t> t;
  ibn::valer<double> E;
  ibn::valer<double> S;
  ibn::valer<double> Ebepc;
  ibn::valer<double> I;
};

template<class T>
inline std::istream & operator >> (std::istream & ifs, ibn::valer<T> & data)
{
  return ifs >> data.value >> data.error;
}

template<class T>
inline std::ostream & operator <<  (std::ostream & ofs, ibn::valer<T> & data)
{
  return ofs << data.value << "+-" <<  data.error;
}


template <class Filter>
std::list<EMS_t>  read_file(std::string file, Filter filter )
{
  std::list<EMS_t>  L;
  std::ifstream ifs(file);
  if(!ifs) { std::cerr << "Unable to open file " << file << std::endl; }
  EMS_t ems;
  while ( ifs >> ems.t >> ems.E >> ems.S >> ems.Ebepc >> ems.I)
  {
    if(!filter(ems)) continue;
    L.push_back(ems);
  }
  return L;
}

TGraphErrors * make_graph(const std::list<EMS_t> & lst)
{
  TGraphErrors *  g = new TGraphErrors;
  for(auto & ems : lst)
  {
    int n = g->GetN();
    g->SetPoint(n, ems.t.value, ems.E.value);
    g->SetPointError(n, ems.t.error, ems.E.error);
  }
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.5);
  return g;
}

int main(int argc, char ** argv)
{
  int nlevels=1;
  double lambda=5000;
  if(argc>=2) nlevels=std::stoi(argv[1]);
  if(argc>=3) lambda = std::stod(argv[2]);

  std::string efile="E.results";
  std::string pfile="P.results";
  std::ifstream ifs(efile);
  if(!ifs) { std::cerr << "Unable to open file " << argv[1] << std::endl; }
  EMS_t ems;
  auto filter = [](const EMS_t & ems){ return ems.t.value>1524e6 && ems.t.error<5000;  };
  auto elst= read_file(efile,filter);
  auto plst= read_file(pfile,filter);
  auto eg = make_graph(elst);
  eg->SetMarkerColor(kBlue);
  eg->SetLineColor(kBlue);

  auto pg = make_graph(plst);
  pg->SetMarkerColor(kRed);
  pg->SetLineColor(kRed);

  TMultiGraph * mg = new TMultiGraph;
  mg->Add(eg, "p");
  mg->Add(pg, "p");
	TApplication theApp("Draw ems files", &argc, argv);
  auto c = new TCanvas;
  c->Divide(1,3);
  c->cd(3);
  mg->Draw("a");
  mg->GetXaxis()->SetTimeDisplay(kTRUE);
  //multistep_fit("efun",eg,nlevels,lambda);
  multistep_fit2("efun","pfun",eg,pg, nlevels,lambda);
  eg->GetFunction("efun")->SetLineColor(kBlack);
  //multistep_fit("pfun",pg, nlevels,lambda);
  pg->GetFunction("pfun")->SetLineColor(kBlack);

  c->cd(1);
  eg->Draw("ap");
  eg->GetXaxis()->SetTimeDisplay(kTRUE);
  c->cd(2);
  pg->Draw("ap");
  pg->GetXaxis()->SetTimeDisplay(kTRUE);
  c->Modified();
  c->Update();

  new TCanvas;
  auto fe = (TF1*)eg->GetFunction("efun")->Clone();
  auto fp = (TF1*)pg->GetFunction("pfun")->Clone();
  std::list<EMS_t> eplst;
  for(const EMS_t & ems: elst) eplst.push_back(ems);
  for(const EMS_t & ems: plst) eplst.push_back(ems);
  eplst.sort([](const auto & e1, const auto & e2) { return e1.t.value < e2.t.value; });
  //std::sort(begin(eplst),end(eplst),  );
  TGraph * fun_g =  new TGraph;
  double tmin = eplst.front().t.value;
  double tmax = eplst.back().t.value;
  for(int i=0;i<1000;i++)
  {
    double t = tmin + (tmax-tmin)*i/1000;
    double E1 = fe->Eval(t);
    double E2 = fp->Eval(t);
    double E = sqrt(E1*E2)*cos(0.011);
    fun_g->SetPoint(fun_g->GetN(), t, E-1776.86);
  }
  //fun_g->SetLineColor(kGreen);
  fun_g->Draw("al");
  fun_g->GetXaxis()->SetTimeDisplay(kTRUE);
  fun_g->GetXaxis()->SetTitle("time");


  theApp.Run();
}



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

inline double multi_step_function(double *x, double * par)
{
  int n = par[0];
  std::vector<double> T(par+1, par+n);
  std::sort(T.begin(),T.end());
  std::vector<double> A(par+n,par+2*n);
  for(int i=0;i<T.size();i++)
  {
    if( x[0] < T[i] ) return A[i];
  }
  return A[T.size()];
}

inline double erf_multi_step_function(double *x, double * par)
{
  int n = par[0];
  double lambda = par[1];
  double xmin = par[2];
  double xmax = par[3];
  std::vector<double> R(par+4, par+4+n); //ranges
  std::vector<double> A(par+n+5,par+2*n+5); //levels
  std::vector<double> T(n-1);
  double range_scale =0;
  for( auto r : R) range_scale += r;
  range_scale = (xmax-xmin)/range_scale;
  double xsum = xmin;
  for(int i=0;i< T.size();i++) 
  {
    xsum += fabs(R[i]*range_scale);
    T[i] = xsum;
  };
  //calculate points
  double a = A[0];
  for(int i=0;i<T.size();i++)
  {
    a+= A[i]*erf((x[0]-T[i])/lambda);
  }
  return a;
}

//inline double erf_multi_step_function2(double *x, double * par)
//{
//  int n = par[0];
//  double lambda = par[1];
//  struct data_t { double t; double a; };
//  std::vector<data_t> D(n);
//  double a = par[2];
//  for(int i=0;i<D.size();i++)
//  {
//    D[i].a = par[3+i];
//    D[i].t = i == (D.size()-1) ?  -std::numeric_limits<double>::infinity() : par[2+n+i];
//  };
//  sort(D.begin(),D.end(), [](const auto &d1, const auto &d2) { return d1.t<d2.t; }
//  //calculate points
//  for(int i=0;i<D.size()-1;i++)
//  {
//    a += (D[i+1].a-D[i].a)*(1.0 + erf((x[0]-D[i])/lambda))*0.5;
//  }
//  return a;
//}


inline void multi_step_fit(const char * name, int nlevels, double reg_lambda, TGraph * g)
{
  TF1 * f = new TF1(name, erf_multi_step_function, 0, 1, 2*nlevels+1);
  int index=0;
  f->FixParameter(index++,nlevels);
  f->FixParameter(index++, reg_lambda);
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  double dx = (xmax-xmin)/nlevels;
  double dy = (ymax-ymin)/nlevels;
  for(int i=0;i<nlevels; i++, index++)  
  {
    f->SetParameter(index, ymin + i*(ymax-ymin)/nlevels);
    f->SetParLimits(index, ymin-dy,ymax+dy);
  }
  for(int i=0;i<nlevels-1; i++, index++) 
  {
    f->SetParLimits(index, xmin,xmax);
    f->SetParameter(index, xmin + (i+1)*dx);
  }
  g->Fit(name);
}

//inline void multi_step_fit3(const char * name, int nlevels, double reg_lambda, TGraph * g)
//{
//  TF1 * f = new TF1(name, erf_multi_step_function2, 0, 1, 2*nlevels+1);
//  int index=0;
//  f->FixParameter(index++,nlevels);
//  f->FixParameter(index++, reg_lambda);
//  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
//  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
//  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
//  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
//  double dx = (xmax-xmin)/nlevels;
//  double dy = (ymax-ymin)/nlevels;
//  for(int i=0;i<nlevels; i++, index++)  
//  {
//    f->SetParameter(index, ymin + i*(ymax-ymin)/nlevels);
//    f->SetParLimits(index, ymin-dy,ymax+dy);
//  }
//  for(int i=0;i<nlevels-1; i++, index++) 
//  {
//    f->SetParLimits(index, xmin,xmax);
//    f->SetParameter(index, xmin + (i+1)*dx);
//  }
//  g->Fit(name);
//}

template <class Filter>
std::list<EMS_t>  read_file(std::string file, Filter filter )
{
  std::list<EMS_t>  L;
  std::ifstream ifs(file);
  if(!ifs) { std::cerr << "Unable to open file " << file << std::endl; }
  EMS_t ems;
  while ( ifs >> ems.t >> ems.E >> ems.S >> ems.Ebepc >> ems.I)
  {
    //std::cout << ems.t << "  " << ems.E << " " << ems.S << " " << ems.Ebepc <<"  " << ems.I << std::endl;
    if(!filter(ems)) continue;
    L.push_back(ems);
  //int i=0;
    //eg->SetPoint(i, ems.t.value, ems.E.value);
    //eg->SetPointError(i, ems.t.error, ems.E.error);
    //i++;
  }
  return L;
}
TGraphErrors * make_graph(const std::list<EMS_t> & lst)
{
  TGraphErrors *  g = new TGraphErrors;
  for(auto & ems : lst)
  {
    int n = g->GetN();
    //std::cout << n << " " << ems.t.value << " " << ems.E.error << std::endl;
    g->SetPoint(n, ems.t.value, ems.E.value);
    g->SetPointError(n, ems.t.error, ems.E.error);
  }
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.5);
  return g;
}

int main(int argc, char ** argv)
{
	TApplication theApp2("Draw ems files", &argc, argv);
  //TF1 * fms = new TF1("fms",multistep2,0,10,3+2*3);
  //fms->SetParameter(0, 3);
  //fms->SetParameter(1, 0.5);
  //fms->SetParameter(2, 0);
  //fms->SetParameter(3, 2);
  //fms->SetParameter(4, 2);
  //fms->SetParameter(5, 5);
  //fms->SetParameter(6, 5);
  //fms->SetParameter(7, 7);
  //fms->SetParameter(8, 7);
  //fms->SetNpx(1000);
  //fms->Draw();
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
  multistep_fit("efun",eg,nlevels,lambda);
  eg->GetFunction("efun")->SetLineColor(kBlack);
  multistep_fit("pfun",pg, nlevels,lambda);
  pg->GetFunction("pfun")->SetLineColor(kBlack);

  c->cd(1);
  eg->Draw("ap");
  c->cd(2);
  pg->Draw("ap");
  c->Modified();
  c->Update();

  //new TCanvas;
  //eg->GetFunction("efun")->Draw();
  //pg->GetFunction("pfun")->Draw("same");
  //double tmin = ltot.front().t.value;
  //double tmax = ltot.back().t.value;
  //TGraph * fun_g =  new TGraph;
  //for(int i=0;i<100;i++)
  //{
  //  double t = tmin + (tmax-tmin)*i/100;
  //  double E1 = eg->GetFunction("efun")->Eval(t);
  //  double E2 = pg->GetFunction("pfun")->Eval(t);
  //  double E = sqrt(E1*E2)*cos(0.011);
  //  fun_g->SetPoint(fun_g->GetN(), t, E);
  //}
  //fun_g->SetLineColor(kGreen);
  //fun_g->Draw("same");
  //theApp2.Run();

  theApp.Run();
}



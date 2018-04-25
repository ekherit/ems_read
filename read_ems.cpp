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

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TF1.h>
#include <TSpline.h>
#include <TSystem.h>


#include <ibn/valer.h>

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

inline double erf_multi_step_function2(double *x, double * par)
{
  int n = par[0];
  double lambda = par[1];
  std::vector<double> T(par+2, par+n+1); //knots
  sort(T.begin(),T.end());
  std::vector<double> A(par+n+1,par+2*n+1); //levels
  //calculate points
  double a = A[0];
  for(int i=0;i<T.size();i++)
  {
    a += (A[i+1]-A[i])*(1.0 + erf((x[0]-T[i])/lambda))*0.5;
  }
  return a;
}

inline void multi_step_fit(const char * name, int nlevels, double reg_lambda, TGraph * g)
{
  TF1 * f = new TF1(name, erf_multi_step_function2, 0, 1, 2*nlevels+1);
  int index=0;
  f->FixParameter(index++,nlevels);
  f->FixParameter(index++, reg_lambda);
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  double dx = (xmax-xmin)/nlevels;
  double dy = (ymax-ymin)/nlevels;
  for(int i=0;i<nlevels-1; i++, index++) 
  {
    f->SetParLimits(index, xmin,xmax);
    f->SetParameter(index, xmin + (i+1)*dx);
  }
  for(int i=0;i<nlevels; i++, index++)  
  {
    f->SetParameter(1+nlevels+i, ymin + i*(ymax-ymin)/(nlevels-1));
    f->SetParLimits(index, ymin-dy,ymax+dy);
  }
  g->Fit(name);
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
    std::cout << ems.t << "  " << ems.E << " " << ems.S << " " << ems.Ebepc <<"  " << ems.I << std::endl;
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
    std::cout << n << " " << ems.t.value << " " << ems.E.error << std::endl;
    g->SetPoint(n, ems.t.value, ems.E.value);
    g->SetPointError(n, ems.t.error, ems.E.error);
  }
  g->SetMarkerStyle(21);
  return g;
}

int main(int argc, char ** argv)
{
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
  multi_step_fit("efun",5,5000,eg);
  eg->GetFunction("efun")->SetLineColor(kBlue);
  multi_step_fit("pfun",5,5000,pg);
  pg->GetFunction("pfun")->SetLineColor(kRed);
	TApplication theApp("Draw ems files", &argc, argv);
  auto c = new TCanvas;
  c->Divide(1,2);

  TGraphErrors * wg = new TGraphErrors;
  TGraphErrors * wg2 = new TGraphErrors;
  std::list<EMS_t> ltot = elst;
  std::list<EMS_t> l2 = plst;
  ltot.merge(l2, [](const auto & p1,const  auto & p2) { return p1.t.value < p2.t.value; });
  ltot.sort([](const auto & p1,const  auto & p2) { return p1.t.value < p2.t.value; });
  

  for(auto & ems :ltot)
  {
    int n = wg->GetN();
    double Ee = eg->Eval(ems.t.value);
    double Ep = pg->Eval(ems.t.value);
    double W = 2.0*sqrt(Ee*Ep)*cos(0.011);
    //std::cout << n << " E = " << Ee << " P = " << Ep << " W=" << W <<std::endl;
    if(std::isnormal(W)) 
    {
      wg->SetPoint(n, ems.t.value, W/2-1776.86);
      wg2->SetPoint(n, ems.t.value, W/2);
    }
  }
  wg2->SetLineWidth(3);
  //mg->Add(wg2,"l");
  c->cd(1);
  mg->Draw("a");
  mg->GetXaxis()->SetTimeDisplay(kTRUE);

  c->cd(2);
  wg->Draw("al");
  wg->GetXaxis()->SetTitle("time");
  wg->GetXaxis()->SetTimeDisplay(kTRUE);
  wg->GetYaxis()->SetTitle("Wcm/2 - M_{#tau}");
  gSystem->ProcessEvents();
  
  std::cout << "Making spline " << std::endl;

  new TCanvas;
  std::vector<double> T;
  std::vector<double> E1;
  std::vector<double> E2;
  double tmin = ltot.front().t.value;
  double tmax = ltot.back().t.value;
  int i=0;
  int Nspline=50;
  while(i<Nspline)
  {
    auto t = tmin + (tmax-tmin)*double(i)/Nspline;
    auto e = wg2->Eval(t);
    if(!std::isnan(e)&& !std::isinf(e)) 
    {
      T.push_back(t);
      E1.push_back(e);
      E2.push_back(e/2-1886.86);
    }
    i++;
  };
  std::cout << "Before splien" << std::endl;
  TSpline3 spline1("s1",T.front(),T.back(),&E1[0],T.size());
  TSpline3 spline2("s2",T.front(),T.back(),&E2[0],T.size());
  c->cd(1);
  spline1.SetLineWidth(2);
  spline1.Draw("same");
  c->cd(2);
  spline2.SetLineColor(kBlue);
  spline2.Draw("same");
  int nlevels=5;
  new TCanvas;
  multi_step_fit("f", 5, 5000,wg);
  theApp.Run();
}



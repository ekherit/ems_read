/*
 * =====================================================================================
 *
 *       Filename:  multistep.h
 *
 *    Description:  Fit data by multistep function
 *
 *        Version:  1.0
 *        Created:  26.04.2018 21:17:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#include <vector>

#include <TF1.h>

inline double multistep(double *x, double * par)
{
  int n = par[0];
  double lambda = par[1];
  // 0 - n
  // 1 - lambda
  // 2 - a0
  // 3 - t0
  // 4 - a1
  // 5 - t1
  // 6 - a2
  // 7 - t2
  // ....
  // 2*n -    an
  // 2*n+1    tn = infinity 
  struct data_t { double a,t; };
  std::vector<data_t> D(n);
  for(int i=0;i<D.size()-1;++i)
  {
    D[i].a = par[2*(i+1)];
    D[i].t = par[2*(i+1)+1];
  }
  D[D.size()-1].a = par[2*n];
  D[D.size()-1].t = std::numeric_limits<double>::max();
  std::sort(D.begin(),D.end(), [](const auto &d1, const auto &d2) { return d1.t<d2.t; });
  //calculate points
  double a = D[0].a;
  for(int i=0;i<D.size()-1;i++)
  {
    a += (D[i+1].a-D[i].a)*(1.0 + erf((x[0]-D[i].t)/lambda))*0.5;
  }
  return a;
}

//multistep with no sorting
//the A is the jumps not the levels except first one(zero one)
inline double multistep2(double *x, double * par)
{
  int n = par[0]; //now n is number of knots - zero meens constant fit 
  double lambda = par[1];
  // 0 - n
  // 1 - lambda
  // 2 - a
  // 3 - t0
  // 4 - a0
  // 5 - t1
  // 6 - a1
  // 7 - t2
  // 8 - a2
  // ....
  // 1+2*n     tn-1
  // 2+2*n -    an-1
  double a = par[2];
  double * knot = &par[3];
  std::vector< std::pair<double, double> > D(n);
  for(int i=0;i<n;++i) 
  {
    D[i].first = knot[2*i]; //knot
    D[i].second = knot[2*i+1]; //level
  }
  std::sort(D.begin(),D.end(), [](const auto &d1, const auto &d2) { return d1.first<d2.first; });
  auto step = [](double xi) { return (1.0 + erf(xi))*0.5; };
  a+=(D[0].second-a)*step((x[0]-D[0].first)/lambda);
  for(int i=1;i<n;i++) a+=(D[i].second-D[i-1].second)*step((x[0]-D[i].first)/lambda);
  return a;
}

inline void multi_step_fit3(const char * name, int nlevels, double reg_lambda, TGraph * g, TCanvas * c)
{
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  double dx = (xmax-xmin)/2;
  double dy = (ymax-ymin)/2;
  //set initial pars
  std::vector<double> pars(2*nlevels+1);
  int index=0;
  pars[index++] = nlevels;
  pars[index++] = reg_lambda;
  //this is for the one knont and two levels
  for(int i=0; i<2; ++i, ++index) pars[index] = ymin + i*(ymax-ymin)/2.;
  pars[index++]=0.5*(xmax+xmin);
  for(; index< pars.size(); ++index) pars[index]=(xmin+xmax)*0.5;
  //now  do the fit
  for(int n = 2; n<=nlevels; n++)
  {
    std::cout << "Doing for " << n << " levels" << std::endl;
    TF1 * f = new TF1(name, multistep, xmin, xmax, 2*n+1);
    f->SetParameters(&pars[0]);
    f->FixParameter(0, n); //fix  current nlevels
    f->FixParameter(1, reg_lambda); //fix lambda
    std::cout << "Setting par lmits" << std::endl;
    for(int i=0;i<n-1;++i) 
    {
      f->SetParLimits(2*(i+1), ymin, ymax);
      f->SetParLimits(2*(i+1)+1, xmin, xmax);
      f->FixParameter(2*(i+1)+1, pars[2*(i+1)+1]);
    }
    f->SetParLimits(2*n, ymin, ymax);

    double chi2 = std::numeric_limits<double>::max();
    double optx;
    std::cout << "Finding optimum" << std::endl;
    for(int i=0;i<g->GetN();++i)
    {
      //std::cout << "Fixing parameters " << 2*n-1 << " to " << g->GetX()[i] << " index=" << i <<std::endl;
      f->FixParameter(2*n-1, g->GetX()[i]);
      g->Fit(name,"Q0");
      double chi2_tmp = g->GetFunction(name)->GetChisquare();
      if( chi2_tmp < chi2)  
      {
        chi2 = chi2_tmp;
        optx = g->GetX()[i];
      }
    }
    std::cout << "Set to found optimum par " << optx << std::endl;
    f->SetParameter(2*n-1,optx);
    std::cout << "Release parameter " << std::endl;
    f->ReleaseParameter(2*n-1);
    g->Fit(name);
    c->Modified();
    c->Update();

    //copy parameters to pars
    f->GetParameters(&pars[0]);
  }
}

inline void multi_step_fit4(const char * name, int nknots, double reg_lambda, TGraph * g, TCanvas * c)
{
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  double dx = (xmax-xmin)/2;
  double dy = (ymax-ymin)/2;
  //set initial pars
  std::vector<double> pars(3+2*nknots);
  int index=0;
  pars[0] = nknots;
  pars[1] = reg_lambda;
  pars[2] = (ymax+ymin)*0.5;
  double * knots  = &pars[3];
  knots[0] = (xmax+xmin)*0.5;
  knots[1] = (ymax-ymin)*0.25;
  for(int n = 1; n<=nknots; ++n)
  {
    std::cout << "Doing for " << n << " levels" << std::endl;
    TF1 * f = new TF1(name, multistep2, xmin, xmax, 3+2*nknots);
    f->SetParameters(&pars[0]);
    f->FixParameter(0, nknots); //fix  current nlevels
    f->FixParameter(1, reg_lambda); //fix lambda
    std::cout << "Setting par lmits" << std::endl;
    for(int i=0;i<n;++i) 
    {
      f->SetParLimits(2*(i+1)+1, xmin, xmax);
      f->SetParLimits(2*(i+1), ymin, ymax);
      f->FixParameter(2*(i+1)+1, pars[2*(i+1)+1]);
    }
    f->SetParLimits(2*n, ymin, ymax);

    double chi2 = std::numeric_limits<double>::max();
    double optx;
    std::cout << "Finding optimum" << std::endl;
    for(int i=0;i<g->GetN();++i)
    {
      //std::cout << "Fixing parameters " << 2*n-1 << " to " << g->GetX()[i] << " index=" << i <<std::endl;
      f->FixParameter(2*n-1, g->GetX()[i]);
      g->Fit(name,"Q0");
      double chi2_tmp = g->GetFunction(name)->GetChisquare();
      if( chi2_tmp < chi2)  
      {
        chi2 = chi2_tmp;
        optx = g->GetX()[i];
      }
    }
    std::cout << "Set to found optimum par " << optx << std::endl;
    f->SetParameter(2*n-1,optx);
    std::cout << "Release parameter " << std::endl;
    f->ReleaseParameter(2*n-1);
    g->Fit(name);
    c->Modified();
    c->Update();

    //copy parameters to pars
    f->GetParameters(&pars[0]);
  }
}

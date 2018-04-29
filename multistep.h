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
  int n         = par[0]; //n - number of levels
  double lambda = par[1];
  // 0 - n
  // 1 - lambda
  // 2 - t0
  // 3 - a0
  // 4 - t1
  // 5 - a1
  // 6 - t2
  // 7 - a2
  // ...
  // 2+2*(n-1)  tn-1
  // 3+2*(n-1)   an-1
  double * knot = &par[2];
  auto step = [](double xi) { return (1.0 + erf(xi))*0.5; };
  double a = knot[0+1];
  for(int i=1;i<n;i++) a+=(knot[2*i+1]-knot[2*(i-1)+1])*step((x[0]-knot[2*i])/lambda);
  return a;
}



double get_chi2(const std::vector<std::pair<double, double> > & knots, int N, double *X, double *Y, double * dX=nullptr, double *dY=nullptr)
{
  double chi2=0;
  for(int i=0;i<N;++i)
  {
    for( int n = 0; n < knots.size()-1; ++n)
    {
      if(knots[n].first < X[i] &&  X[i] < knots[n+1].first)
      {
        double error = 1;
        if( dY != nullptr )
        {
          if( dY[i]!= 0 ) error = dY[i];
        }
        chi2+=pow((knots[n].second - Y[i])/error,2.0);
      }
    }
    //special case last knot
    if(knots[knots.size()-1].first < X[i])
    {
        double error = 1;
        if( dY != nullptr )
        {
          if( dY[i]!= 0 ) error = dY[i];
        }
        chi2+=pow((knots[knots.size()-1].second - Y[i])/error,2.0);
    }
  }
  return chi2;
}

double fit_knots(std::vector<std::pair<double, double> > & knots, int N, double *X, double *Y, double * dX = nullptr, double * dY = nullptr)
{
  std::sort(begin(knots), end(knots), [](const auto & p1, const auto &p2){ return p1.first<p2.first;} );
  std::vector<long> count(knots.size());
  for(auto & p : knots) p.second = 0;
  for(auto & c: count) c=0;

  for(int i=0;i<N;++i)
  {
    for( int n = 0; n < knots.size()-1; ++n)
    {
      if(knots[n].first < X[i] &&  X[i] < knots[n+1].first)
      {
        count[n]++;
        knots[n].second += Y[i];
      }
    }
    //special case last knot
    if(knots[knots.size()-1].first < X[i])
    {
        count[knots.size()-1]++;
        knots[knots.size()-1].second += Y[i];
    }
  }
  for(int n=0;n<knots.size();++n)
  {
    if(count[n]!=0) knots[n].second/=count[n];
  }
  return get_chi2(knots, N, X, Y, dX, dY);
}

std::vector<std::pair<double, double> > estimate_knots(int nlevels, int N, double * X, double *Y, double *dX = nullptr, double *dY = nullptr)
{
  std::vector<std::pair<double, double> > knots;
  knots.reserve(nlevels);
  knots.push_back({-std::numeric_limits<double>::max(), Y[0]});
  for(int n=1;n<nlevels; ++n )
  {
    knots.push_back({X[0],Y[N-1]});
    double x=X[0];
    double chi2=std::numeric_limits<double>::max();
    for( int i=0;i<N;++i)
    {
      knots[n].first = X[i];
      std::vector<std::pair<double,double>> tmp_knots = knots;
      double tmp_chi2 = fit_knots(tmp_knots, N,X,Y,dX,dY);
      if(tmp_chi2<chi2)
      {
        chi2 = tmp_chi2;
        x=X[i];
      }
    }
    knots[n].first = x;
    fit_knots(knots, N,X,Y,dX,dY);
  }
  return knots;
}

typedef std::vector<std::pair<double, double> > knots_t;
std::vector<knots_t>  estimate_knots2(
    int nlevels, 
    int N1, double * X1, double *Y1,
    int N2, double * X2, double *Y2
    )
{
  std::vector<std::pair<double, double> > knots1,knots2;
  knots1.reserve(nlevels);
  knots2.reserve(nlevels);
  knots1.push_back({-std::numeric_limits<double>::max(), Y1[0]});
  knots2.push_back({-std::numeric_limits<double>::max(), Y2[0]});
  for(int n=1;n<nlevels; ++n )
  {
    knots1.push_back({std::min(X1[0],X2[0]),Y1[N1-1]});
    knots2.push_back({std::min(X1[0],X2[0]),Y2[N2-1]});
    double chi2=std::numeric_limits<double>::max();

    double x=X1[0];
    for( int i=0;i<N1;++i)
    {
      knots1[n].first = X1[i];
      knots2[n].first = X1[i];
      std::vector<std::pair<double,double>> tmp_knots1 = knots1;
      std::vector<std::pair<double,double>> tmp_knots2 = knots2;
      double tmp_chi2=0;
      tmp_chi2 += fit_knots(tmp_knots1, N1,X1,Y1);
      tmp_chi2 += fit_knots(tmp_knots2, N2,X2,Y2);
      if(tmp_chi2<chi2)
      {
        chi2 = tmp_chi2;
        x=X1[i];
      }
    }
    for( int i=0;i<N2;++i)
    {
      knots1[n].first = X2[i];
      knots2[n].first = X2[i];
      std::vector<std::pair<double,double>> tmp_knots1 = knots1;
      std::vector<std::pair<double,double>> tmp_knots2 = knots2;
      double tmp_chi2=0;
      tmp_chi2 += fit_knots(tmp_knots1, N1,X1,Y1);
      tmp_chi2 += fit_knots(tmp_knots2, N2,X2,Y2);
      if(tmp_chi2<chi2)
      {
        chi2 = tmp_chi2;
        x=X2[i];
      }
    }
    knots1[n].first = x;
    knots2[n].first = x;
    fit_knots(knots1, N1,X1,Y1);
    fit_knots(knots2, N2,X2,Y2);
  }
  return {knots1, knots2};
}

TF1 * make_knot_function(const char * name, const std::vector< std::pair<double, double> > & knots, double lambda, double xmin, double xmax, double ymin, double ymax)
{
  TF1 * f = new TF1(name, multistep, 0, 1, 2 + 2*knots.size());
  f->FixParameter(0, knots.size());
  f->FixParameter(1, lambda);
  for( int i=0;i< knots.size();++i)
  {
    f->SetParameter(2+2*i, knots[i].first);
    f->SetParameter(2+2*i+1, knots[i].second);
  }
  f->FixParameter(2,knots[0].first);
  for( int i=1;i<knots.size();++i) f->SetParLimits(2+2*i, xmin, xmax);
  for( int i=0;i<knots.size();++i) f->SetParLimits(2+2*i+1, ymin, ymax);
  f->SetNpx(1000);
  return f;
}

TF1 * multistep_fit(const char * name, TGraph * g, double nlevels, double lambda=1)
{
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto kn = estimate_knots(nlevels, g->GetN(), g->GetX(), g->GetY());
  TF1 * f = make_knot_function(name, kn, lambda, xmin, xmax, ymin, ymax);
  g->Fit(name);
  return f;
}

std::vector<TF1*> multistep_fit2(const char * name1, const char *name2, TGraph * g1, TGraph * g2, double nlevels, double lambda=1)
{
  auto xmin1 = *std::min_element(&g1->GetX()[0], &g1->GetX()[g1->GetN()-1]);
  auto xmax1 = *std::max_element(&g1->GetX()[0], &g1->GetX()[g1->GetN()-1]);
  auto ymin1 = *std::min_element(&g1->GetY()[0], &g1->GetY()[g1->GetN()-1]);
  auto ymax1 = *std::max_element(&g1->GetY()[0], &g1->GetY()[g1->GetN()-1]);

  auto xmin2 = *std::min_element(&g2->GetX()[0], &g2->GetX()[g2->GetN()-1]);
  auto xmax2 = *std::max_element(&g2->GetX()[0], &g2->GetX()[g2->GetN()-1]);
  auto ymin2 = *std::min_element(&g2->GetY()[0], &g2->GetY()[g2->GetN()-1]);
  auto ymax2 = *std::max_element(&g2->GetY()[0], &g2->GetY()[g2->GetN()-1]);
  auto xmin = std::min(xmin1, xmin2);
  auto xmax = std::max(xmax1, xmax2);
  auto ymin = std::min(ymin1, ymin2);
  auto ymax = std::max(ymax1, ymax2);

  auto kn = estimate_knots2(nlevels, 
      g1->GetN(), g1->GetX(), g1->GetY(),
      g2->GetN(), g2->GetX(), g2->GetY()
      );
  TF1 * f1 = make_knot_function(name1, kn[0], lambda, xmin, xmax, ymin, ymax);
  TF1 * f2 = make_knot_function(name2, kn[1], lambda, xmin, xmax, ymin, ymax);
  for(int i=0;i<nlevels;i++) 
  {
    f1->FixParameter(2+2*i, kn[0][i].first);
    f2->FixParameter(2+2*i, kn[0][i].first);
  }
  g1->Fit(name1);
  g2->Fit(name2);
  return {f1,f2};
}

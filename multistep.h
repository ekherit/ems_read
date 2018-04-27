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
  int n         = par[0]; //now n is number of knots - zero meens constant fit
  double lambda = par[1];
  double a      = par[2];
  // 0 - n
  // 1 - lambda
  // 2 - 
  // 3 - t0
  // 4 - a0
  // 5 - t1
  // 6 - a1
  // 7 - t2
  // 8 - a2
  // ....
  // 1+2*n     tn-1
  // 2+2*n -    an-1
  double * knot = &par[3];
  std::vector< std::pair<double, double> > D(n);
  //std::cout << "Initial level : " << a << std::endl;
  for(int i=0;i<n;++i) 
  {
    D[i].first = knot[2*i]; //knot
    D[i].second = knot[2*i+1]; //level
    //std::cout << " knot " << i << " " << D[i].first << "   level " << D[i].second << std::endl;
  }
  std::sort(D.begin(),D.end(), [](const auto &d1, const auto &d2) { return d1.first<d2.first; });
  auto step = [](double xi) { return (1.0 + erf(xi))*0.5; };
  a+=(D[0].second-a)*step((x[0]-D[0].first)/lambda);
  for(int i=1;i<n;i++) a+=(D[i].second-D[i-1].second)*step((x[0]-D[i].first)/lambda);
  return a;
}

inline double multistep3(double *x, double * par)
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
  int Npars = 3+2*nknots;
  std::vector<double> pars(Npars);
  //std::vector<double> knots(nknots);
  int index=0;
  pars[0] = nknots;
  pars[1] = reg_lambda;
  pars[2] = (ymax+ymin)*0.5;
  double * knots  = &pars[3];
  knots[0] = (xmax+xmin)*0.5;
  knots[1] = (ymax-ymin)*0.25;
  std::cout << "Before knot cycle " << std::endl;
  for(int n = 1; n<=nknots; ++n)
  {
    std::cout << "Doing for " << n << " knots" << std::endl;
    TF1 * f = new TF1(name, multistep2, xmin, xmax, 3+2*nknots);
    //f->SetParameters(&pars[0]); //copy parameters from previouse
    //initital function paramete set
    f->FixParameter(0, nknots); //fix  current nlevels
    f->FixParameter(1, reg_lambda); //fix lambda
    f->SetParLimits(2, ymin, ymax);
    for(int i=0;i<n;++i)
    {
      f->SetParameter(2*(i+1)+1, 0.5*(xmax+xmin));
      f->SetParameter(2*(i+1)+2, 0.5*(ymax+ymin));
      f->SetParLimits(2*(i+1)+1, xmin, xmax);
      f->SetParLimits(2*(i+1)+2, xmin, xmax);
    }
    double chi2 = std::numeric_limits<double>::max();
    double optx;
    //knots[2*(n-1)] = x;
    std::vector< std::pair<double, double> > D(n);
    int x_index = 0;
    for(int i=0;i<g->GetN();++i)
    {
      double x = g->GetX()[i];
      std::cout << "n=" << n << " i=" << i << "/" << g->GetN() <<  "  " << (x-g->GetX()[0])/60. << " npars = " << Npars  << std::endl;
      pars[3+2*n] = x; //add knots to result pars
      for(int i=0;i<n;++i) 
      {
        D[i].first = knots[2*i]; //knot
        D[i].second = knots[2*i+1]; //level
        //std::cout << " knot " << i << " " << D[i].first << "   level " << D[i].second << std::endl;
      }
      //sort parameters
      std::sort(D.begin(),D.end(), [](const auto &d1, const auto &d2) { return d1.first<d2.first; });
      for(int i=0;i<n;++i) std::cout << (D[i].first-g->GetX()[0])/60. << " ";
      std::cout << std::endl;
      //fix and set all new parameters layout
      //int left_free_index, right_free_index;
      int tmp_idx=0;
      for(int j=0;j<n;++j) 
      {
        //knots[2*i] = D[i].first;
        //knots[2*i+1] = D[i].second;
        f->FixParameter(3+2*j,   D[j].first);
        f->SetParameter(3+2*j+1, D[j].second);
        //f->FixParameter(3+2*j+1, D[j].second);
        //if( x > D[j].first  && j<n-1 ) 
        //{
        //  f->ReleaseParameter(3+2*j+1);
        //  f->ReleaseParameter(3+2*(j+1)+1);
        //}
        if( x == D[j].first )  tmp_idx = j;
      }
      //g->Fit(name,"Q0");
      g->Fit(name);
      double chi2_tmp = g->GetFunction(name)->GetChisquare();
      if( chi2_tmp < chi2 )  
      {
        chi2 = chi2_tmp;
        x_index = tmp_idx;
        optx = g->GetX()[i];
      }
    }
    std::cout << "x_index=" << x_index << " optx = " << optx << std::endl;
    f->SetParameter(3+2*x_index, optx);
    f->ReleaseParameter(3+2*x_index);
    g->Fit(name);
    //c->Modified();
    //c->Update();
    //copy parameters to pars
    f->GetParameters(&pars[0]);
    std::cout << "Before next knot" << std::endl;
  }
}

double get_chi2(const std::vector<std::pair<double, double> > & knots, int N, double *X, double *Y)
{
  //sort knots firtly
  //auto sorted_knots = knots;
  //std::sort(begin(sorted_knots), end(sorted_knots), [](const auto & p1, const auto &p2){ return p1.first<p2.first;} );
  double chi2=0;
  for(int i=0;i<N;++i)
  {
    for( int n = 0; n < knots.size()-1; ++n)
    {
      if(knots[n].first < X[i] &&  X[i] < knots[n+1].first)
      {
        chi2+=pow(knots[n].second - Y[i],2.0);
      }
    }
    //special case last knot
    if(knots[knots.size()-1].first < X[i])
    {
        chi2+=pow(knots[knots.size()-1].second - Y[i],2.0);
    }
  }
  return chi2;
}

double fit_knots(std::vector<std::pair<double, double> > & knots, int N, double *X, double *Y)
{
  //for(auto & p : knots)
  //{
  //std::cout << "fit_knots before sort: " << (p.first-X[0])/60.0 << "  " << p.second << std::endl;
  //}
  std::sort(begin(knots), end(knots), [](const auto & p1, const auto &p2){ return p1.first<p2.first;} );
  std::vector<long> count(knots.size());
  //for(auto & p : knots)
  //{
  //std::cout << "fit_knots: " << (p.first-X[0])/60.0 << "  " << p.second << std::endl;
  //}
  for(auto & p : knots) p.second = 0;
  for(auto & c: count) c=0;

  for(int i=0;i<N;++i)
  {
    //std::cout << "examine " << i << (X[i]-X[0])/60.0 << "  " << Y[i] << std::endl;
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
  //for(int i=0;i<count.size(); i++) std::cout << i << " " << count[i] << std::endl;
  return get_chi2(knots, N, X, Y);
}


std::vector<std::pair<double, double> > estimate_knots(int nlevels, int N, double * X, double *Y)
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
      double tmp_chi2 = fit_knots(tmp_knots, N,X,Y);
      if(tmp_chi2<chi2)
      {
        chi2 = tmp_chi2;
        x=X[i];
      }
    }
    knots[n].first = x;
    fit_knots(knots, N,X,Y);
  }
  return knots;
}

TF1 * make_knot_function(const char * name, const std::vector< std::pair<double, double> > & knots, double lambda, double xmin, double xmax)
{
  TF1 * f = new TF1(name, multistep3, 0, 1, 2 + 2*knots.size());
  f->FixParameter(0, knots.size());
  f->FixParameter(1, lambda);
  for( int i=0;i< knots.size();++i)
  {
    f->SetParameter(2+2*i, knots[i].first);
    f->SetParameter(2+2*i+1, knots[i].second);
  }
  f->FixParameter(2,knots[0].first);
  for( int i=1;i<knots.size();++i) f->SetParLimits(2+2*i, xmin, xmax);
  return f;
}

TF1 * multistep_fit(const char * name, TGraph * g, double nlevels, double lambda=1)
{
  auto xmin = *std::min_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto xmax = *std::max_element(&g->GetX()[0], &g->GetX()[g->GetN()-1]);
  auto ymin = *std::min_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto ymax = *std::max_element(&g->GetY()[0], &g->GetY()[g->GetN()-1]);
  auto kn = estimate_knots(nlevels, g->GetN(), g->GetX(), g->GetY());
  TF1 * f = make_knot_function(name, kn, lambda, xmin, xmax);
  g->Fit(name);
  return f;
}

#include "interp.h"

MatrixXd AtA(MatrixXd A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
    .rankUpdate(A.adjoint());
}

double threshold(){
  return 1.0E6; // ???? FIXME
}

// see https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf, Fig. 9:
ArrayXd Dplus(const ArrayXd& d) {
  ArrayXd di(d.size());
  double comp(d.maxCoeff() * threshold());
  for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
  return di;
}

double kern2d(double x, double xi, double hx,
		     double y, double yi, double hy,
		     std::string kernel){
  // implement product kernels
  double t1, t2, k;

  if(kernel=="gaussian"){
    // hx is interpreted as 3*sx ... so
    hx=hx/3.0;
    hy=hy/3.0;
  }

  t1=(x-xi)/hx;
  t2=(y-yi)/hy;

  //Rcout << "t1: " << t1 << " t2: " << t2;
  if(kernel=="gaussian")
    k=1.0/(2.0*M_PI)*exp(-0.5*(t1*t1+t2*t2));
  else if(kernel=="epanechnikov"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0))
      k=3.0*3.0/4.0/4.0*(1-t1*t1)*(1-t2*t2);
    else
      k=0.0;
  } else if(kernel=="biweight"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0))
      k=15.0*15.0/16.0/16.0*(1-t1*t1)*(1-t1*t1)*(1-t2*t2)*(1-t2*t2);
    else
      k=0.0;
  } else if(kernel=="tricube"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0)){
      double t1a=abs(t1), t2a=abs(t2);
      k=70.0*70.0/81.0/81.0*(1-t1a*t1a*t1a)*(1-t1a*t1a*t1a)*(1-t1a*t1a*t1a)*(1-t2a*t2a*t2a)*(1-t2a*t2a*t2a)*(1-t2a*t2a*t2a);
    } else
      k=0.0;
  } else if(kernel=="triweight"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0))
      k=35.0*35.0/32.0/32.0*(1-t1*t1)*(1-t1*t1)*(1-t1*t1)*(1-t2*t2)*(1-t2*t2)*(1-t2*t2);
    else
      k=0.0;
  } else if(kernel=="cosine"){
    if((abs(t1)<=M_PI/2.0) && (abs(t2)<=M_PI/2.0))
      k=0.25*cos(t1)*cos(t2);
    else
      k=0.0;
  } else if(kernel=="uniform"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0))
      k=0.25;
    else
      k=0.0;
  } else if(kernel=="triangle"){
    if((abs(t1)<=1.0) && (abs(t2)<=1.0))
      k=(1.0-abs(t1))*(1.0-abs(t2));
    else
      k=0.0;
  } else
    Rf_error("kernel not implemented!");
  //Rcout << " k: " << k << std::endl;
  return k;
}

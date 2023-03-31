#include <Rcpp.h>
// fuer string Vergleiche:
#include <string>

using namespace Rcpp;

#include "interp.h"

/* Implementation according to
 * [1] Akima, H. (1970) A new method of interpolation
 *     and smooth curve fitting based on local procedures,
 *     J. ACM \bold{17}(4), 589-602
 *
 * [2] Akima, H. (1991) A Method of Univariate Interpolation that Has
 *     the Accuracy of a Third-degree Polynomial. ACM Transactions on
 *     Mathematical Software, \bold{17}(3), 341-366.
 */

// [2] eqn (8) whith notation xji=(xj-xi), xki=(xk-xi), xli=(xl-xi)
#define F0(xji,xki,xli,yji,yki,yli)             \
  (  yji*xki*xki*xli*xli*(xli-xki)              \
  + yki*xli*xli*xji*xji*(xji-xli)               \
  + yli*xji*xji*xki*xki*(xki-xji)               \
  / (xji*xki*xli*(xki-xji)*(xli-xki)*(xli-xji))

#define F(i,j,k,l) ( (y[j]-y[i]) * (x[k]-x[i])*(x[k]-x[i]) * (x[l]-x[i])*(x[l]-x[i]) * (x[l]-x[k]) \
                     + (y[k]-y[i]) * (x[l]-x[i])*(x[l]-x[i]) * (x[j]-x[i])*(x[j]-x[i]) * (x[j]-x[l]) \
                     + (y[l]-y[i]) * (x[j]-x[i])*(x[j]-x[i]) * (x[k]-x[i])*(x[k]-x[i]) * (x[k]-x[j]) ) \
  / ( (x[j]-x[i])*(x[k]-x[i])*(x[l]-x[i])*(x[k]-x[j])*(x[l]-x[k])*(x[l]-x[j]) )

// [2] eqn (11)
#define D(i,j,k,l) (x[j]-x[i])*(x[j]-x[i]) + (x[k]-x[i])*(x[k]-x[i]) + (x[l]-x[i])*(x[l]-x[i])

// [2] sum y^2
#define S(i,j,k,l) (y[i]*y[i] + y[j]*y[j] + y[k]*y[k] + y[l]*y[l])

// [2] eqn (10)
#define b0(i,j,k,l) ((x[i]*x[i]+x[j]*x[j]+x[k]*x[k]+x[l]*x[l])*(y[i]+y[j]+y[k]+y[l])-(x[i]+x[j]+x[k]+x[l])*(x[i]*y[i]+x[j]*y[j]+x[k]*y[k]+x[l]*y[l])) \
  / (4.0*(x[i]*x[i]+x[j]*x[j]+x[k]*x[k]+x[l]*x[l]) - (x[i]+x[j]+x[k]+x[l])*(x[i]+x[j]+x[k]+x[l]))

#define b1(i,j,k,l) (4.0*(x[i]*y[i]+x[j]*y[j]+x[k]*y[k]+x[l]*y[l])-(x[i]+x[j]+x[k]+x[l])*(y[i]+y[j]+y[k]+y[l])) \
  / (4.0*(x[i]*x[i]+x[j]*x[j]+x[k]*x[k]+x[l]*x[l]) - (x[i]+x[j]+x[k]+x[l])*(x[i]+x[j]+x[k]+x[l]))

// [2] eqn (9)
#define V(i,j,k,l) (y[i]-(b0(i,j,k,l)+b1(i,j,k,l)*x[i]))*(y[i]-(b0(i,j,k,l)+b1(i,j,k,l)*x[i])) \
  + (y[j]-(b0(i,j,k,l)+b1(i,j,k,l)*x[j]))*(y[j]-(b0(i,j,k,l)+b1(i,j,k,l)*x[j])) \
  + (y[k]-(b0(i,j,k,l)+b1(i,j,k,l)*x[k]))*(y[k]-(b0(i,j,k,l)+b1(i,j,k,l)*x[k])) \
  + (y[l]-(b0(i,j,k,l)+b1(i,j,k,l)*x[l]))*(y[l]-(b0(i,j,k,l)+b1(i,j,k,l)*x[l]))

// [[Rcpp::export()]]
List aSpline(NumericVector x, 
             NumericVector y,
             NumericVector xout,
             CharacterVector method="improved",
             int degree=3
             ) {

  List ret;
  
  int nx = x.size();
  NumericVector yp=NumericVector(nx);
  // check for xout=NULL done in R wrapper!
  int n=xout.size();    
  
  
  NumericVector yout=NumericVector(n);

  double a0, a1, a2, a3, mj;
  double ypmmj, ypmj, yppj, ypppj;
  double wpmmj, wpmj, wppj, wpppj;
  double x2, x3, y2, y3;
  double eps=1.0e-12,f,v,d,xmin,xmax;

  xmin=x[0];
  xmax=x[nx-1];
  
  for(int j=0; j<nx; j++){
    double xm1=R_NegInf;
    // check for increasing data, get min/max x
    if(j==0)
      xm1=x[j];
    else{
      if(x[j]<xm1)
        Rf_error("points in data set not in increasing order!");
      if(x[j]==xm1)
        Rf_error("points in data set repeated!");
      xm1=x[j];        
    }
    if(x[j]<xmin)
      xmin=x[j];
    if(x[j]>xmax)
      xmax=x[j];
  }


  // estimate partial derivatives for method [2] 
  if(as<std::string>(method)=="improved"){
    for(int j=0; j<nx; j++){
   
      if(j>2 && j<nx){
        f = F(j,j-3,j-2,j-1);
        ypmmj = f;
        v=V(j,j-3,j-2,j-1);
        d=D(j,j-3,j-2,j-1);
        //Rcout <<"ijkl: " << j << j-3 << j-2 << j-1 << " F: " << f<< " V: " << v << " D: " << d << std::endl;
        if(v>eps*S(j,j-3,j-2,j-1))
          wpmmj=1.0/(v*d);
        else
          wpmmj=1.0;
      } else {
        ypmmj=0;
        wpmmj=0;
      }
      //Rcout <<"ijkl: " << j << j-3 << j-2 << j-1 << " wpmmj: " << wpmmj << std::endl;
    
      if(j>1 && j<nx-2){
        f = F(j,j-2,j-1,j+1);
        ypmj = f;
        v=V(j,j-2,j-1,j+1);
        d=D(j,j-2,j-1,j+1);
        //Rcout <<"ijkl: " << j << j-2 << j-1 << j+1<< " F: " << f << " V: " << v << " D: " << d << std::endl;
        if(v>eps*S(j,j-2,j-1,j+1))
          wpmj=1.0/(v*d);
        else
          wpmj=0;
      } else {
        ypmj = 0;
        wpmj=0;
      }
      //Rcout <<"ijkl: " << j << j-2 << j-1 << j+1<< " ypmj: " << ypmj << " wpmj: " << wpmj << std::endl;

      if(j>0 && j<nx-3){
        f = F(j,j-1,j+1,j+2);
        yppj = f;
        v=V(j,j-1,j+1,j+2);
        d=D(j,j-1,j+1,j+2);
        //Rcout <<"ijkl: " << j << j-1 << j+1 << j+2<< " F: " << f << " V: " << v << " D: " << d << std::endl;
        if(v>eps*S(j,j-1,j+1,j+2))
          wppj=1.0/(v*d);
        else
          wppj=1.0;
      } else {
        yppj=0;
        wppj=0;
      }
      //Rcout <<"ijkl: " << j << j-1 << j+1 << j+2<< " yppj: " << yppj << " wppj: " << wppj << std::endl;
      if(j<nx-4){
        f= F(j,j+1,j+2,j+3);
        ypppj= f;
        v=V(j,j+1,j+2,j+3);
        d=D(j,j+1,j+2,j+3);
        //Rcout <<"ijkl: " << j << j+1 << j+2 << j+3<< " F: " << f << " V: " << v << " D: " << d << std::endl;
        if(v>eps*S(j,j+1,j+2,j+3))
          wpppj=1.0/(v*d);
        else
          wpppj=1.0;
      } else {
        ypppj=0;
        wpppj=0;
      }
      //Rcout <<"ijkl: " << j << j+1 << j+2 << j+3<< " wpppj: " << wpppj << std::endl;
      // build average estimate of first derivative:
      yp[j]=(ypmmj*wpmmj+ypmj*wpmj+yppj*wppj+ypppj*wpppj)/(wpmmj+wpmj+wppj+wpppj);
      //Rcout << "final yp: " << yp[j] << std::endl;
    }
  }

  
  // for method [1] add two estimated points on a parabola below and above data:
  // from [1] eqn (8) and [1] eqn (9)
  double xm1=0.0,xm2=0.0,xnp0=0.0,xnp1=0.0, // means: x[-1], x[-2], x[n+0], x[n+1]
    ym1=0.0,ym2=0.0,ynp0=0.0,ynp1=0.0, //        y[-1], y[-2], y[n+0], y[n+1]
    g0l=0.0,g1l=0.0,g2l=0.0,g0u=0.0,g1u=0.0,g2u=0.0; // coefficients of parabolas for extra points below/above range of x
  
  if(as<std::string>(method)=="original"){
    /*
    if(nx>5){
      Rcout << "x: " << x[0] << " " << x[1] << " " << x[2] << " ... " << x[nx-3] << " " << x[nx-2] << " " << x[nx-1] << std::endl;
      Rcout << "y: " << y[0] << " " << y[1] << " " << y[2] << " ... " << y[nx-3] << " " << y[nx-2] << " " << y[nx-1] << std::endl;
    } else if(nx==5){
      Rcout << "x: " << x[0] << " " << x[1] << " " << x[2] << " " << x[nx-2] << " " << x[nx-1] << std::endl;
      Rcout << "y: " << y[0] << " " << y[1] << " " << y[2] << " " << y[nx-2] << " " << y[nx-1] << std::endl;
    } else if(nx==4){
      Rcout << "x: " << x[0] << " " << x[1] << " " << x[nx-2] << " " << x[nx-1] << std::endl;
      Rcout << "y: " << y[0] << " " << y[1] << " " << y[nx-2] << " " << y[nx-1] << std::endl;
    } else if(nx==3){
      Rcout << "x: " << x[0] << " " << x[1] << " " << x[2] << std::endl;
      Rcout << "y: " << y[0] << " " << y[1] << " " << y[2]  << std::endl;
    } else if(nx==2){
      Rcout << "x: " << x[0] << " " << x[1] << std::endl;
      Rcout << "y: " << y[0] << " " << y[1] << std::endl;
    }else if(nx==1){
      Rcout << "x: " << x[0] << " " << std::endl;
      Rcout << "y: " << y[0] << " " << std::endl;
    }
    */
    // extra points below/above x range:
    // nx=2 needs extra handling as there is no x[2]
    if(nx==2){
      xm1  = x[0]   -x[1]   +x[0];
      xm2  = xm1   - x[0]   +xm1;
      xnp0 = x[nx-1]-x[nx-2]+x[nx-1];
      xnp1 = xnp0   -x[nx-1]+xnp0;
      ym1 =y[0]   - (x[0]-xm1)    *((y[1]-y[0])/(x[1]-x[0]));
      ym2 =y[0]   - (x[0]-xm2)    *((y[1]-y[0])/(x[1]-x[0]));
      ynp0=y[nx-1]+ (xnp0-x[nx-1])*((y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2]));
      ynp1=y[nx-1]+ (xnp1-x[nx-1])*(2.0*(ynp0-y[nx-1])/(xnp0-x[nx-1]));
    } else {
      xm1  = x[1]   -x[2]   +x[0];
      xm2  = x[0]   -x[1]   +xm1;
      xnp0 = x[nx-1]-x[nx-3]+x[nx-2];
      xnp1 = xnp0   -x[nx-2]+x[nx-1];
      ym1 =y[0]   + (x[0]-xm1)    *((y[2]-y[1])/(x[2]-x[1])-2.0*(y[1]-y[0])/(x[1]-x[0]));
      ym2 =ym1    + (xm1-xm2)     *((y[1]-y[0])/(x[1]-x[0])-2.0*(y[0]-ym1)/(x[0]-xm1));
      ynp0=y[nx-1]+ (xnp0-x[nx-1])*(2.0*(y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2])-(y[nx-2]-y[nx-3])/(x[nx-2]-x[nx-3]));
      ynp1=ynp0   + (xnp1-xnp0)   *(2.0*(ynp0-y[nx-1])/(xnp0-x[nx-1])-(y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2]));
      // save parameters of g_0+g_1 x+g_2 x^2 for later reuse during extrapolation,
      // needed twice:
      //
      if(nx>2){
	g0l = (x[0]*x[1]*x[1]*y[2]-x[0]*x[0]*x[1]*y[2]-x[0]*x[2]*x[2]*y[1]+x[0]*x[0]*x[2]*y[1]
	       +x[1]*x[2]*x[2]*y[0]-x[1]*x[1]*x[2]*y[0])/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]));
	g1l = -(x[1]*x[1]*y[2]-x[0]*x[0]*y[2]-x[2]*x[2]*y[1]+x[0]*x[0]*y[1]+x[2]*x[2]*y[0]-x[1]*x[1]*y[0])/
	  ((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]));
	g2l = (x[1]*y[2]-x[0]*y[2]-x[2]*y[1]+x[0]*y[1]+x[2]*y[0]-x[1]*y[0])/
	  ((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]));
	g0u = (x[nx-3]*x[nx-2]*x[nx-2]*y[nx-1]-x[nx-3]*x[nx-3]*x[nx-2]*y[nx-1]-x[nx-3]*x[nx-1]*x[nx-1]*y[nx-2]+x[nx-3]*x[nx-3]*x[nx-1]*y[nx-2]
	       +x[nx-2]*x[nx-1]*x[nx-1]*y[nx-3]-x[nx-2]*x[nx-2]*x[nx-1]*y[nx-3])/((x[nx-2]-x[nx-3])*(x[nx-1]-x[nx-3])*(x[nx-1]-x[nx-2]));
	g1u = -(x[nx-2]*x[nx-2]*y[nx-1]-x[nx-3]*x[nx-3]*y[nx-1]-x[nx-1]*x[nx-1]*y[nx-2]+x[nx-3]*x[nx-3]*y[nx-2]+x[nx-1]*x[nx-1]*y[nx-3]-x[nx-2]*x[nx-2]*y[nx-3])/
	  ((x[nx-2]-x[nx-3])*(x[nx-1]-x[nx-3])*(x[nx-1]-x[nx-2]));
	g2u = (x[nx-2]*y[nx-1]-x[nx-3]*y[nx-1]-x[nx-1]*y[nx-2]+x[nx-3]*y[nx-2]+x[nx-1]*y[nx-3]-x[nx-2]*y[nx-3])/
	  ((x[nx-2]-x[nx-3])*(x[nx-1]-x[nx-3])*(x[nx-1]-x[nx-2]));
      } else if(nx==2){
	g0l = -(x[0]*y[1]-x[1]*y[0])/(x[1]-x[0]);
	g1l = (y[1]-y[0])/(x[1]-x[0]);
	g2l = 0.0;
	g0u = g0l;
	g1u = g1l;
	g2u = 0.0;
	//Rcout << "g0l:" << g0l << " g1l:" << g1l << " g2l:" << g2l << std::endl;
	//Rcout << "g0u:" << g0u << " g1u:" << g1u << " g2u:" << g2u << std::endl;
      }
    }
    
    /*
    Rcout << "size x: " << x.size() << std::endl;
    Rcout << "size y: " << y.size() << std::endl;
    if(nx==2){
      Rcout << "extra x: " << xm2 << ", " << xm1 << ", (" << x[0] << ")...(" << x[nx-1] << "), " << xnp0 << ", " << xnp1 << std::endl;
      Rcout << "extra y: " << ym2 << ", " << ym1 << ", (" << y[0] << ")...(" << y[nx-1] << "), " << ynp0 << ", " << ynp1 << std::endl;
    } else if(nx==3){
      Rcout << "extra x: " << xm2 << ", " << xm1 << ", (" << x[0] << "),("<<x[1]<<")...(" << x[nx-1] << "), " << xnp0 << ", " << xnp1 << std::endl;
      Rcout << "extra y: " << ym2 << ", " << ym1 << ", (" << y[0] << "),("<<y[1]<<")...(" << y[nx-1] << "), " << ynp0 << ", " << ynp1 << std::endl;
    } else if(nx==4){
      Rcout << "extra x: " << xm2 << ", " << xm1 << ", (" << x[0] << "),("<<x[1]<<"),("<<x[2]<<")...(" << x[nx-1] << "), " << xnp0 << ", " << xnp1 << std::endl;
      Rcout << "extra y: " << ym2 << ", " << ym1 << ", (" << y[0] << "),("<<y[1]<<"),("<<y[2]<<")...(" << y[nx-1] << "), " << ynp0 << ", " << ynp1 << std::endl;
    } else {
      Rcout << "extra x: " << xm2 << ", " << xm1 << ", (" << x[0] << "),("<<x[1]<<"),("<<x[2]<<"),("<<x[3]<<")...(" << x[nx-1] << "), " << xnp0 << ", " << xnp1 << std::endl;
      Rcout << "extra y: " << ym2 << ", " << ym1 << ", (" << y[0] << "),("<<y[1]<<"),("<<y[2]<<"),("<<y[3]<<")...(" << y[nx-1] << "), " << ynp0 << ", " << ynp1 << std::endl;
    }
    */

  }


  // iterate over output points
  for(int i=0; i<n; i++){
    // special cases n=1,2,3,4, mostly for method [2] because
    // method [1] uses the extra points below and above xmin and xmax
    if(nx==1){
      // stupid case: single data point -> constant output
      yout[i]=y[0];
      Rf_warning("only one point in data set!");
    } else if(nx==2){// && as<std::string>(method)=="improved"){
      // two data points -> straight line
      yout[i]=y[0]+(y[1]-y[0])/(x[1]-x[0])*(xout[i]-x[0]);
    } else if(nx==3 && as<std::string>(method)=="improved"){
      // three points: quadratic function with linear extrapolation
      
      // transform x[i] -> x[i]-x[0], y[i] -> y[0]
      // => a0=0 and x1=0, y1=0

      //x1=x[0]-x[0]=0;
      x2=x[1]-x[0];
      x3=x[2]-x[0];
      //y1=y[0]-y[0]=0;
      y2=y[1]-y[0];
      y3=y[2]-y[0];
      double denom=x2*(x2-x3)*x3;
      if(std::fabs(denom)<1.0e-16){
        Rf_error("points in data set coincide!");
      }
      // solution to a0+a1*xi+a2*xi^2=yi, i=1,2,3
      // a0=0.0;
      a1=(-x3*x3*y2+x2*x2*y3)/denom;
      a2=(x3*y2-x2*y3)/denom;
      //Rcout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << std::endl;
      // transform back
      // extrapolation with linear function
      if(xout[i]<xmin){
        yout[i]=y[0]+ a1 *(xout[i]-x[0]);
      } else if(xout[i]>xmax){
        yout[i]=y[2]+ (a1+2.0*a2*x3) *(xout[i]-x[2]);
      } else {
        // interpolation
        yout[i]=y[0]+(xout[i]-x[0])*(a1+(xout[i]-x[0])*a2);
      }
    } else if(nx==4 && as<std::string>(method)=="improved"){
      // four points: cubic function with linear extrapolation
      
      // transform x[i] -> x[i]-x[0], y[i] -> y[0]
      // => a0=0 and x1=0, y1=0
      
      double //x1=x[0]-x[0]=0,
        x2=x[1]-x[0],
        x3=x[2]-x[0],
        x4=x[3]-x[0],
        //y1=y[0]-y[0]=0,
        y2=y[1]-y[0],
        y3=y[2]-y[0],
        y4=y[3]-y[0];
      double denom=x2*(x2-x3)*x3*(x2-x4)*(x3-x4)*x4;
      if(std::fabs(denom)<eps){
        Rf_error("points in data set coincide!");
      }
      // solution to a0+a1*xi+a2*xi^2+a3*xi^3=yi, i=1,2,3,4
      // a0=0.0;
      a1=(x2*x2*x4*x4*(x4-x2)*y3
          +x3*x3*x3*(x4*x4*y2-x2*x2*y4)
          +x3*x3*(x2*x2*x2*y4-x4*x4*x4*y2)
          )/denom;
      a2=(x2*x4*y3*(x2*x2-x4*x4)
          +x3*x3*x3*(x2*y4-x4*y2)
          +x3*(x4*x4*x4*y2-x2*x2*x2*y4)
          )/denom;
      a3=(x2*x4*y3*(x4-x2)
          +x3*x3*(x4*y2-x2*y4)
          +x3*(x2*x2*y4-x4*x4*y2)
          )/denom;

      //Rcout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << " a3: " << a3 << std::endl;
      // transform back
      // extrapolation (linear function)
      if(xout[i]<xmin){
        yout[i]=y[0]+ a1 *(xout[i]-x[0]);
      } else if(xout[i]>xmax){
        yout[i]=y[3]+ (a1+2.0*a2*x4+3.0*a3*x4*x4) *(xout[i]-x[3]);
      } else {
        // interpolation
        yout[i]=y[0]+(xout[i]-x[0])*(a1+(xout[i]-x[0])*(a2+(xout[i]-x[0])*a3));
      }
    } else { // main case: method="original" or more then 4 points:
      if(as<std::string>(method)=="improved"){
	// extrapolation
	if(xout[i]<xmin){
	  yout[i]=y[0]+yp[0]*(xout[i]-x[0]);
        } else if(xout[i]>=xmax){
	  yout[i]=y[nx-1]+yp[nx-1]*(xout[i]-x[nx-1]);
	} else {
	  // interpolation ([2])
          for(int j=0; j<nx-1; j++){
            if(x[j]<=xout[i] && xout[i]<x[j+1]){
              // [2] eqn (5):
              mj=(y[j+1]-y[j])/(x[j+1]-x[j]);
              // [2] eqn (4)
              a0=y[j];
              a1=yp[j]; 
              a2=-(2.0*(yp[j]-mj)+yp[j+1]-mj)/(x[j+1]-x[j]);
              a3=(yp[j]-mj+yp[j+1]-mj)/((x[j+1]-x[j])*(x[j+1]-x[j]));
              //Rcout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << " a3: " << a3 << std::endl;
              if(degree==3){
                // numerically instable:
                // yout[i]=a0+a1*(xout[i]-x[j])+a1*(xout[i]-x[j])*(xout[i]-x[j])+a2*(xout[i]-x[j])*(xout[i]-x[j])*(xout[i]-x[j]);
                // better:
                yout[i]=a0+(xout[i]-x[j])*(a1+(xout[i]-x[j])*(a2+(xout[i]-x[j])*a3));
              } else {
                // TODO a0*(u^d-u)+a1*((1-u)^d-(1-u))
                // transformation x,y -> u,v [2] eqn (15)
                double ui=(xout[i]-x[j])/(x[j+1]-x[j]);
                
                // [2] eqn (17):
                mj=(y[j+1]-y[j])/(x[j+1]-x[j]);
                //Rcout << "mj: " << mj << " ypj: " << yp[j] << " ypj+1: " << yp[j+1] << std::endl;
                // [2] eqn [18]
                double v0p=(yp[j]-mj)*(x[j+1]-x[j]),
                  v1p=(yp[j+1]-mj)*(x[j+1]-x[j]);
                // v(u)=A0(u^n-u)+A1((1-u)^n-(1-u))
                // [2] eqn (20)
                double A0=(v0p+((double)degree-1)*v1p)/((double)degree*((double)degree-2));
                double A1=-(((double)degree-1)*v0p+v1p)/((double)degree*((double)degree-2));
		
                //Rcout << "A: " << A0 << " A1: " << A1 << std::endl;
                double vi=A0*(std::pow(ui,degree)-ui)+A1*(std::pow(1-ui,degree)-(1-ui));
                //Rcout << "ui: " << ui << " vi: " << vi << std::endl;
		
                yout[i]=vi+y[j]+(y[j+1]-y[j])*ui;
                //Rcout << "yout: " << yout[i] << std::endl;
              }
            }
          }
	}
      } else if(as<std::string>(method)=="original"){
	double m1=0.0, m2=0.0, m3=0.0, m4=0.0, m5=0.0, t1, t2;
	double x1=0.0, x2=0.0, y1=0.0, y2=0.0, p0, p1, p2, p3;
	// extrapolation
	if(xout[i]<xmin){
	  //Rcout << "extrapolation <xmin" << std::endl;
          //x1=x[0]; x2=x[1]; y1=y[0]; y2=y[1];
	  yout[i]=g0l+(xout[i])*(g1l+(xout[i])*g2l);
	  /*
          // slopes from five points, according to p590 in [1]          
          m1=(ym1-ym2)/(xm1-xm2); 
          m2=(y[0]-ym1)/(x[0]-xm1); 
          m3=(y[1]-y[0])/(x[1]-x[0]);
	  if(nx>=4){
	  m4=(y[2]-y[1])/(x[2]-x[1]);  
	  m5=(y[3]-y[2])/(x[3]-x[2]);
	  } else if(nx==3){
	  m4=(y[2]-y[1])/(x[2]-x[1]);  
	  m5=(ynp0-y[2])/(xnp0-x[2]);
	  } else if(nx==2){
	  m4=(ynp0-y[1])/(xnp0-x[1]);  
	  m5=(ynp1-ynp0)/(xnp1-xnp0);
	  }
	  */
	} else if(xout[i]>=xmax){
	  //Rcout << "extrapolation >xmax" << std::endl;
          //x1=x[nx-2]; x2=x[nx-1]; y1=y[nx-2]; y2=y[nx-1];
	  yout[i]=g0u+(xout[i])*(g1u+(xout[i])*g2u);
	  /*
          // slopes from five points, according to p590 in [1]
	  if(nx>=4){
	  m1=(y[nx-3]-y[nx-4])/(x[nx-3]-x[nx-4]); 
	  m2=(y[nx-2]-y[nx-3])/(x[nx-2]-x[nx-3]); 
	  } else if(nx==3){
	  m1=(y[nx-3]-ym1)/(x[nx-3]-xm1); 
	  m2=(y[nx-2]-y[nx-3])/(x[nx-2]-x[nx-3]); 
	  } else if(nx==2){
	  m1=(ym1-ym2)/(xm1-xm2); 
	  m2=(y[nx-2]-ym1)/(x[nx-2]-xm1); 
	  }
          m3=(y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2]);  
          m4=(ynp0-y[nx-1])/(xnp0-x[nx-1]);  
          m5=(ynp1-ynp0)/(xnp1-xnp0); 
	  */       
	} else {
	  //Rcout << "interpolate" << std::endl;
          
	  for(int j=0; j<nx-1; j++){
	    if(x[j]<=xout[i] && xout[i]<x[j+1]){

	      x1=x[j], x2=x[j+1], y1=y[j], y2=y[j+1];
	      // slopes from five points, according to p590 in [1]
	      if(j==0){
		m1=(ym1-ym2)/(xm1-xm2); 
		m2=(y[0]-ym1)/(x[0]-xm1); 
		m3=(y[1]-y[0])/(x[1]-x[0]);  
		m4=(y[2]-y[1])/(x[2]-x[1]);  
		m5=(y[3]-y[2])/(x[3]-x[2]); 
	      } else if(j==1){
		m1=(y[0]-ym1)/(x[0]-xm1); 
		m2=(y[1]-y[0])/(x[1]-x[0]); 
		m3=(y[2]-y[1])/(x[2]-x[1]);  
		m4=(y[3]-y[2])/(x[3]-x[2]);  
		m5=(y[4]-y[3])/(x[4]-x[3]); 
	      } else if(j+1==nx-2){
		m1=(y[nx-4]-y[nx-5])/(x[nx-4]-x[nx-5]); 
		m2=(y[nx-3]-y[nx-4])/(x[nx-3]-x[nx-4]); 
		m3=(y[nx-2]-y[nx-3])/(x[nx-2]-x[nx-3]); 
		m4=(y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2]);  
		m5=(ynp0-y[nx-1])/(xnp0-x[nx-1]);  
	      } else if(j+1==nx-1){
		m1=(y[nx-3]-y[nx-4])/(x[nx-3]-x[nx-4]); 
		m2=(y[nx-2]-y[nx-3])/(x[nx-2]-x[nx-3]); 
		m3=(y[nx-1]-y[nx-2])/(x[nx-1]-x[nx-2]);  
		m4=(ynp0-y[nx-1])/(xnp0-x[nx-1]);  
		m5=(ynp1-ynp0)/(xnp1-xnp0); 
	      } else {
		m1=(y[j-1]-y[j-2])/(x[j-1]-x[j-2]); 
		m2=(y[j]-y[j-1])/(x[j]-x[j-1]); 
		m3=(y[j+1]-y[j])/(x[j+1]-x[j]);  
		m4=(y[j+2]-y[j+1])/(x[j+2]-x[j+1]);  
		m5=(y[j+3]-y[j+2])/(x[j+3]-x[j+2]); 
	      }
            }
	  }
	  //Rcout << "ms: " << m1 << " " << m2<< " "  << m3<< " "  << m4<< " "  << m5 << std::endl;
	  // [1] eqn (1)
	  // [1] eqn (1)
          if(std::fabs(m1-m2)<eps && std::fabs(m3-m4)>eps)
            t1=m2;
          else if(std::fabs(m1-m2)>eps && std::fabs(m3-m4)<eps)
            t1=m3;
          else if(std::fabs(m1-m2)<eps && std::fabs(m3-m4)<eps)
            t1=0.5*(m2+m3);
          else
            t1=(std::fabs(m4-m3)*m2+std::fabs(m2-m1)*m3)/(std::fabs(m4-m3)+std::fabs(m2-m1));
          
          if(std::fabs(m2-m3)<eps && std::fabs(m4-m5)>eps)
            t2=m3;
          else if(std::fabs(m2-m3)>eps && std::fabs(m4-m5)<eps)
            t2=m4;
          else if(std::fabs(m2-m3)<eps && std::fabs(m4-m5)<eps)
            t2=0.5*(m3+m4);
          else
            t2=(std::fabs(m5-m4)*m3+std::fabs(m3-m2)*m4)/(std::fabs(m5-m4)+std::fabs(m3-m2));

	    
	  //Rcout << "t1: " << t1 << " t2: " << t2 << std::endl;
	  // [1] eqn (3) to (6)
	  p0=y1;
	  p1=t1;
	  p2=(3.0*(y2-y1)/(x2-x1)-2.0*t1-t2)/(x2-x1);
	  p3=(t1+t2-2.0*(y2-y1)/(x2-x1))/((x2-x1)*(x2-x1));
	  //Rcout << "p0:" << p0 << " p1:" << p1 << " p2:" << p2 << " p3:" << p3 << std::endl;
	  // [1] eqn [2]
	  yout[i]=p0+(xout[i]-x1)*(p1+(xout[i]-x1)*(p2+(xout[i]-x1)*p3));
	}

	/*	    
		    if(xout[i]<xmin){
		    yout[i]=g0l+(xout[i])*(g1l+(xout[i])*g2l);
		    Rcout << "g0l:" << g0l << " g1l:" << g1l << " g2l:" << g2l << std::endl;
		    } else if(xout[i]>xmax){
		    yout[i]=g0u+(xout[i])*(g1u+(xout[i])*g2u);
		    Rcout << "g0u:" << g0u << " g1u:" << g1u << " g2u:" << g2u << std::endl;
		    }
	*/
	//Rcout << "out (x,y) = (" << xout[i] << " ," << yout[i] << ")" << std::endl;
      } else {
	Rf_error("unknown method!");
      }
      
    }
    //Rcout << "i: " << i << " xout[i]: " << xout[i] << " yput[i]: " << yout[i] << std::endl;
  }
  ret=List::create(_("x")=xout, _("y")=yout);
  return ret;
}


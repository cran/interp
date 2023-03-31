#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
List BiLinear(NumericVector x, NumericVector y, NumericMatrix z, NumericVector x0, NumericVector y0) {
  List ret;
  int nx=x.size();
  int ny=y.size();
  int n0=x0.size();
  NumericVector z0=NumericVector(n0); 
  if(n0!=y0.size()){
    Rf_error("sizes of x0 and y0 differ!");
  }

  for(int k=0;k<n0;k++)
    for(int i=0;i<nx-1;i++)
      for(int j=0;j<ny-1;j++){
	double x1,y1,xt,yt;
	if(x(i)<=x0(k) and x0(k)<=x(i+1) and
	   y(j)<=y0(k) and y0(k)<=y(j+1)){
 	  x1=x(i+1)-x(i);
	  y1=y(j+1)-y(j);
	  if(x1==0.0 or y1==0.0){
	    Rf_error("some grid step size is zero!");
	  }
	  xt=(x0(k)-x(i))/x1;
	  yt=(y0(k)-y(j))/y1;
	  z0(k)=(1.0-yt)*(1.0-xt)*z(i,j)+
	    (1.0-yt)*xt*z(i+1,j)+
	    yt*(1.0-xt)*z(i,j+1)+
	    yt*xt*z(i+1,j+1); 
	}
      }
  ret=List::create(_("x0")=x0, _("y0")=y0, _("z0")=z0);
  return ret;
}


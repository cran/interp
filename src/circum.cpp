

#include "interp.h"


// [[Rcpp::export]]
List circum(NumericVector x, NumericVector y){

  int nx = x.size();
  int ny = y.size();
  List ret;

  if(nx!=ny)
    Rf_error("size of x and y differs!");

  try {
    double a,b,c;
    
    if(x[0]==x[1] && y[0]==y[1]){
      Rf_error("point 1 and 2 coincide!");
    } else {
      a=sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
    }
    if(x[1]==x[2] && y[1]==y[2]){
      Rf_error("point 2 and 3 coincide!");
    } else {
      b=sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    }
    if(x[2]==x[0] && y[2]==y[0]){
      Rf_error("point 3 and 1 coincide!");
    } else {
      c=sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
    }

    float sp, area, cr, ir, ar;
    // semiperimeter
    sp=(a+b+c)/2.0;
    // area (Heron):
    area=sqrt(sp*(sp-a)*(sp-b)*(sp-c));
    // circumcircle radius
    cr=a*b*c/(4.0*area);
    // inscribed circle radius
    ir=area/sp;
    // aspect ratio
    ar=ir/cr;
    // barycentric coordinates of circumcircle center
    double p,q,r,s;
    p=a*a*(-a*a+b*b+c*c);
    q=b*b*(a*a-b*b+c*c);
    r=c*c*(a*a+b*b-c*c);
    s=p+q+r;
    p=p/s; q=q/s;r=r/s;
    // p associated with side a opposite to x[2],y[2]
    // q       ...            b   ...       x[0],y[0]
    // r       ...            c   ...       x[1],y[1]
    // cartesian coordinates of circumcircle center
    double xc, yc;
    xc=q*x[0]+r*x[1]+p*x[2];
    yc=q*y[0]+r*y[1]+p*y[2];
    // orientation, signed area, positive if counter clockwise
    double sa;
    int orient;
    orient=sgn((y[1]-y[0])*(x[2]-x[1])-(x[1]-x[0])*(y[2]-y[1]));
    sa=orient*area;
    ret=List::create(_("x")=xc, _("y")=yc,
		     _("aspect.ratio")=ar,
		     _("x")=xc,
		     _("y")=yc,
		     _("radius")=cr,
		     _("signed.area")=sa);

    return ret;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached
}
  

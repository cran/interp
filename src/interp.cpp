

#include "interp.h"


// [[Rcpp::export]]
List interpDeltri(NumericVector x, NumericVector y,
                     NumericVector zD,
                     List t, // data xD and yD contained here!
                     CharacterVector input = "points",
                     CharacterVector output = "grid") {

  List T(t);
  int nT = T.size();
  int nG = x.size();
  int mG = y.size();

  NumericMatrix z;
  // initialize return matrix with NA:
  if(as<std::string>(output)=="grid"){
    NumericMatrix z = NumericMatrix(nG,mG,NumericVector (nG*mG,NumericVector::get_na()).begin());
  }
  if(as<std::string>(output)=="points"){
    NumericMatrix z = NumericMatrix(nG,1,NumericVector (nG,NumericVector::get_na()).begin());
  }

  List ret;

  // bounding box for triangles:
  IntegerVector jTsw(nT);
  IntegerVector kTsw(nT);
  IntegerVector jTne(nT);
  IntegerVector kTne(nT);

  try {

    if(as<std::string>(output)=="grid"){
      // get bounding boxes (SW <-> NE) for all triangles:
      for(int i=0; i<nT; i++) {
        SEXP Ti = T[i];

        DataFrame Triangle(Ti);

        NumericVector xT = Triangle["x"];
        NumericVector yT = Triangle["y"];
        NumericVector zT = Triangle["z"];

        // bounding box for triangle i
        double xsw=min(xT);
        double ysw=min(yT);
        double xne=max(xT);
        double yne=max(yT);

        // translate bounding box into grid indices
        jTsw[i]=0;
        kTsw[i]=0;
        jTne[i]=nG-1;
        kTne[i]=mG-1;

        for(int j=0; j<nG; j++){
          if(x[j]<xsw) jTsw[i]=j;
          if(x[nG-j-1]>xne) jTne[i]=nG-j-1;
        }
        for(int k=0; k<mG; k++){
          if(y[k]<ysw) kTsw[i]=1;
          if(y[mG-k-1]>yne) kTne[i]=mG-k-1;
        }
      }
    }

    // iterate over triangles
    for(int i=0; i<nT; i++) {
      SEXP Ti = T[i];

      DataFrame Triangle(Ti);

      NumericVector xT = Triangle["x"];
      NumericVector yT = Triangle["y"];
      NumericVector zT = Triangle["z"];

      if(as<std::string>(output)=="grid"){
        // iterate only over grid points (j,k) inside bounding box of triangle i
        for(int j=jTsw[i]; j<jTne[i]; j++) {
          for(int k=kTsw[i]; k<kTne[i]; k++) {
            // calculate barycentric coordinates:
            double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[k] - yT[2])) /
              ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
            double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[k] - yT[2])) /
              ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
            double c = 1 - a - b;
            // check if inside triangle, handle only yet untouched grid points
            if(R_IsNA(z(j,k))){
              if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
                // perform barycentric interpolation:
                z(j,k)=a*zT[0]+b*zT[1]+c*zT[2];
              }
            }
          }
        }
      } else if(as<std::string>(output)=="points"){
        // iterate over output points
        for(int j=0; j<nG; j++) {
          // calculate barycentric coordinates:
          double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[j] - yT[2])) /
            ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
          double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[j] - yT[2])) /
            ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
          double c = 1 - a - b;
          // check if inside triangle, handle only yet untouched grid points
          if(R_IsNA(z(j,0))){
            if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
              // perform barycentric interpolation:
              z(j,0)=a*zT[0]+b*zT[1]+c*zT[2];
            }
          }
        }
      } else Rf_error("invalid output specification!");
    }

    ret=List::create(_("x")=x, _("y")=y, _("z")=z);

    return ret ;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached

}


// [[Rcpp::export]]
List interpShull(NumericVector x, NumericVector y,
                 NumericVector xD, NumericVector yD,
                 NumericVector zD,
                 bool linear=true,
                 CharacterVector input = "points",
                 CharacterVector output = "grid",
                 CharacterVector kernel = "gaussian",
                 NumericVector h = NumericVector::create(0.0),
                 CharacterVector solver = "QR",
                 int degree = 3,
                 bool baryweight = true,
                 bool autodegree = false,
                 double adtol = 1E-6,
                 bool smoothpde = false,
                 bool akimaweight = true,
                 int nweight = 25) {

  int nxD=xD.size();
  int nyD=yD.size();

  if(as<std::string>(input)=="points" && nxD!=nyD)
    ::Rf_error("length of xD and yD dont match!");

  int nG = x.size();
  int mG = y.size();


  NumericMatrix z;

  // initialize return matrix with NA:
  if(as<std::string>(output)=="grid"){
    z = NumericMatrix(nG,mG,NumericVector (nG*mG,NumericVector::get_na()).begin());
  }
  if(as<std::string>(output)=="points"){
    z = NumericMatrix(nG,1,NumericVector (nG,NumericVector::get_na()).begin());
  }

  List ret;

  try{

    // part 0

    // do s-Hull triangulation:
    // call shDt
    Triang tXY=shDt(Rcpp::as<std::vector<double> >(xD),
		    Rcpp::as<std::vector<double> >(yD));
    // note: triangles are enumerated counterclockwise



    int nT=tXY.nT;

    // Rcout << "get bounding boxes" <<std::endl;

    // get bounding boxes (SW <-> NE) for all triangles:

    IntegerVector jTsw(nT);
    IntegerVector kTsw(nT);
    IntegerVector jTne(nT);
    IntegerVector kTne(nT);

    IntegerVector iT = IntegerVector(3);
    NumericVector xT = NumericVector(3);
    NumericVector yT = NumericVector(3);
    NumericVector zT = NumericVector(3);

    for(int i=0; i<nT; i++) {

      iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
      xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
      yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];



      // bounding box for triangle i
      double xsw=min(xT);
      double ysw=min(yT);
      double xne=max(xT);
      double yne=max(yT);

      // translate bounding box into grid indices, start with complete grid:
      jTsw[i]=0;
      kTsw[i]=0;
      jTne[i]=nG-1;
      kTne[i]=mG-1;

      // Rcout << "bb to grid indices ..." << std::endl;

      for(int j=0; j<nG; j++){
	if(x[j]<xsw) jTsw[i]=j;
	if(x[nG-j-1]>xne) jTne[i]=nG-j-1;
      }
      for(int k=0; k<mG; k++){
	if(y[k]<ysw) kTsw[i]=1;
	if(y[nG-k-1]>yne) kTne[i]=mG-k-1;
      }
    }

    // part 1
    // TODO for !linear

    // Prepare edge structure for discontinuity checks:
    // nCheck points per edge, per triangle side, EdgeCheck.nE is number of edges
    int nCheck=3;
    Edges EdgeCheck;
    MatrixXi eij(nxD,nxD);
    if(!linear){
      EdgeCheck.i1  = VectorXi(nT*nCheck);  // nT*nCheck is surely > nE
      EdgeCheck.i2  = VectorXi(nT*nCheck);
      EdgeCheck.t1  = VectorXi(nT*nCheck);
      EdgeCheck.t2  = VectorXi(nT*nCheck);
      EdgeCheck.xB  = MatrixXd(nT*nCheck,3);
      EdgeCheck.yB  = MatrixXd(nT*nCheck,3);
      EdgeCheck.zBl = MatrixXd(nT*nCheck,3);
      EdgeCheck.zBr = MatrixXd(nT*nCheck,3);

      for(int e=0; e<nT*nCheck; e++){
        EdgeCheck.i1[e]=-2;
        EdgeCheck.i2[e]=-2;
        EdgeCheck.t1[e]=-2;
        EdgeCheck.t2[e]=-2;
      }

      int e=0;
      for(int i=0; i<nxD-1; i++){
        for(int j=i+1; j<nxD; j++){
          EdgeCheck.i1[e]=-1;
          EdgeCheck.i2[e]=-1;
          EdgeCheck.t1[e]=-1;
          EdgeCheck.t2[e]=-1;

          for(int k=0; k<nT; k++){
            if((xD[i]==xD[tXY.i1[k]] && yD[i]==yD[tXY.i1[k]] &&
                ((xD[j]==xD[tXY.i2[k]] && yD[j]==yD[tXY.i2[k]]) ||
                 (xD[j]==xD[tXY.i3[k]] && yD[j]==yD[tXY.i3[k]]))) ||
               (xD[i]==xD[tXY.i2[k]] && yD[i]==yD[tXY.i2[k]] &&
                ((xD[j]==xD[tXY.i1[k]] && yD[j]==yD[tXY.i1[k]]) ||
                 (xD[j]==xD[tXY.i2[k]] && yD[j]==yD[tXY.i2[k]]))) ||
               (xD[i]==xD[tXY.i3[k]] && yD[i]==yD[tXY.i3[k]] &&
                ((xD[j]==xD[tXY.i1[k]] && yD[j]==yD[tXY.i1[k]]) ||
                 (xD[j]==xD[tXY.i2[k]] && yD[j]==yD[tXY.i2[k]]))))
              {
                // take nCheck points between (xD[i],yD[i]) and (xD[j],yD[j])

                if(EdgeCheck.t1[e]==-1){
                  for(int l=0;l<nCheck;l++){
                    double gamma=(1.0*(l+1))/(nCheck+1);
                    EdgeCheck.xB(e,l)=gamma*xD[i]+(1-gamma)*xD[j];
                    EdgeCheck.yB(e,l)=gamma*yD[i]+(1-gamma)*yD[j];
                  }

                  EdgeCheck.t1[e]=k;
                  EdgeCheck.i1[e]=i;
                  EdgeCheck.i2[e]=j;
                  eij(i,j)=e;
                  e++;
                } else
                  EdgeCheck.t2[eij(i,j)]=k;
              }
          }
        }
      }
      EdgeCheck.nE=e;
    }

    // indicator vector fo estimating / re-estimating:
    IntegerVector doEstD=IntegerVector(nxD); // per data point
    NumericMatrix doEstT=NumericMatrix(nT,3); // per triangle side midpoints
    
    // start with local polynomials with given degree
    if(!linear){
      for(int i=0; i<nxD; i++)
        doEstD[i]=degree;
    }

    // prepare variables for estimating partial derivatives
    NumericVector Z(nxD), Z_x(nxD), Z_y(nxD), Z_xy(nxD), Z_xx(nxD), Z_yy(nxD),
      Z_x_ab(nT),Z_y_ab(nT),Z_x_bc(nT),Z_y_bc(nT),Z_x_ca(nT),Z_y_ca(nT);

    int p=0;
    PDEst pde;
    NN nn;
    // Rcout << "start interpolation in " << nT << " triangles" << std::endl;
    VectorXd par_o(21);
    MatrixXd par(nT*3,21); // spline parameters
    // spline parameters:
    double a00,a01,a02,a03,a04,a05,
      a10,a11,a12,a13,a14,
      a20,a21,a22,a23,
      a30,a31,a32,
      a40,a41,
      a50;

    MatrixXd AfTr(nT*3,9);  // to keep all affine transformation matrices, for all 3
    // enumerations
    MatrixXd A(3,3);  // a single affine transformation matrix

    
    if(!linear){

      if(degree==0)
        p=1; // local constant trend
      else if(degree==1)
        p=3; // local linear trend
      else if(degree==2)
        p=6; // local quadratic trend
      else if(degree==3)
        p=10; // local cubic trend
      else if(degree>3)
        Rf_error("degree>3 !");

      // get local neigbouhood structure
      nn=nN(xD,yD);

      for(int dg=degree;dg>=0;dg--){
        for(int i=0; i<nxD; i++){
          if(doEstD[i]<=dg){
	    if(smoothpde)
	      pde=pDsmooth(xD,yD,zD,nn,xD[i],yD[i],kernel,h,as<std::string>(solver),doEstD[i],nweight,akimaweight);
	    else
	      pde=pD(xD,yD,zD,nn,xD[i],yD[i],kernel,h,as<std::string>(solver),doEstD[i]);

            // double sse2=pde.se.transpose()*pde.se;

            Z[i]  =  pde.est[0];
            if(dg>=1){
              Z_x[i]  = pde.est[1];
              Z_y[i]  = pde.est[2];
            } else {
              Z_x[i]  = 0.0;
              Z_y[i]  = 0.0;
            }
            if(dg>=2){
              Z_xy[i] =    pde.est[3];
              Z_xx[i] =    pde.est[4];
              Z_yy[i] =    pde.est[5];
            } else {
              Z_xy[i] = 0.0;
              Z_xx[i] = 0.0;
              Z_yy[i] = 0.0;
            }

          } // if estimate
          //sumDoEstD += doEstD[i];
        } // loop over data
        //if(!autodegree) // leave while loop:
        //  sumDoEstD=-nxD;


        // iterate over triangles
        for(int i=0; i<nT; i++) {


          // ???? Average over all three representations of the triangle,
          // as always (0,0) has errors almost 0
          // as step one, store results for all representations:
          int barycycles=1;
          if(baryweight)
            barycycles=3;

          for(int o=0;o<barycycles;o++){
            if(o==0){
              iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
              xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
              yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];
              zT[0]=zD[tXY.i1[i]]; zT[1]=zD[tXY.i2[i]]; zT[2]=zD[tXY.i3[i]];
            } else if(o==1){
              iT[0]=tXY.i2[i]; iT[1]=tXY.i3[i]; iT[2]=tXY.i1[i];
              xT[0]=xD[tXY.i2[i]]; xT[1]=xD[tXY.i3[i]]; xT[2]=xD[tXY.i1[i]];
              yT[0]=yD[tXY.i2[i]]; yT[1]=yD[tXY.i3[i]]; yT[2]=yD[tXY.i1[i]];
              zT[0]=zD[tXY.i2[i]]; zT[1]=zD[tXY.i3[i]]; zT[2]=zD[tXY.i1[i]];
            } else {
              iT[0]=tXY.i3[i]; iT[1]=tXY.i1[i]; iT[2]=tXY.i2[i];
              xT[0]=xD[tXY.i3[i]]; xT[1]=xD[tXY.i1[i]]; xT[2]=xD[tXY.i2[i]];
              yT[0]=yD[tXY.i3[i]]; yT[1]=yD[tXY.i1[i]]; yT[2]=yD[tXY.i2[i]];
              zT[0]=zD[tXY.i3[i]]; zT[1]=zD[tXY.i1[i]]; zT[2]=zD[tXY.i2[i]];
            }




            double xab,xbc,xca,yab,ybc,yca;

 
            // eqns for the gradient normal to ab, bc, ca:
            // xab=1/2*(xa+xb), yab=1/2*(ya+yb) ...
            xab=0.5*(xT[0]+xT[1]);
            yab=0.5*(yT[0]+yT[1]);
            pde=pD(xD,yD,zD,nn,xab,yab,kernel,h,as<std::string>(solver),degree); // Akima uses 3rd order polynom for local fit

            // TODO: include into reestimate cycle
            if(p>1){
              Z_x_ab[i]  = pde.est[1]; // Note: this is stored twice
              //                         (in the neighbour triangle)!
              Z_y_ab[i]  = pde.est[2];
            } else {
              Z_x_ab[i]  = 0.0;
              Z_y_ab[i]  = 0.0;
            }
            xbc=0.5*(xT[1]+xT[2]);
            ybc=0.5*(yT[1]+yT[2]);
            pde=pD(xD,yD,zD,nn,xbc,ybc,kernel,h,as<std::string>(solver),degree);

            if(p>1){
              Z_x_bc[i]  = pde.est[1];
              Z_y_bc[i]  = pde.est[2];
            } else {
              Z_x_bc[i]  = 0.0;
              Z_y_bc[i]  = 0.0;
            }

            xca=0.5*(xT[2]+xT[0]);
            yca=0.5*(yT[2]+yT[0]);
            pde=pD(xD,yD,zD,nn,xca,yca,kernel,h,as<std::string>(solver),degree);
            /*
              if((pde.se(9)/pde.se(0))>100 || (pde.se(8)/pde.se(0))>100)
              pde=pD(xD,yD,zD,nn,xca,yca,kernel,h,solver,2); // fall back to 2nd order
            */
            if(p>1){
              Z_x_ca[i]  = pde.est[1];
              Z_y_ca[i]  = pde.est[2];
            } else {
              Z_x_ca[i]  = 0.0;
              Z_y_ca[i]  = 0.0;
            }



            // local copies of xa,yz ...

            double xa,ya,xb,yb,xc,yc;

            xa=xT(0);
            xb=xT(1);
            xc=xT(2);
            ya=yT(0);
            yb=yT(1);
            yc=yT(2);


            // parameters for affine transformation:
            /* in Maxima:

               solve A x = t, x=(xa,ya,1), (xb,yb,1), (xc,yc,1) t=(0,0), (0,1), (1,0)
               for aij



               A:matrix([a11,a12,a13],[a21,a22,a23],[a31,a32,a33]);

               in homogenous coordniates:

               a:[xa,ya,1];
               b:[xb,yb,1];
               c:[xc,yc,1];
               ta:[0,0,1];
               tb:[1,0,1];
               tc:[0,1,1];

               matrix equations,row by row:

               ea1:(A.a)[1][1]=ta[1];
               ea2:(A.a)[2][1]=ta[2];
               ea3:(A.a)[3][1]=ta[3];

               eb1:(A.b)[1][1]=tb[1];
               eb2:(A.b)[2][1]=tb[2];
               eb3:(A.b)[3][1]=tb[3];

               ec1:(A.c)[1][1]=tc[1];
               ec2:(A.c)[2][1]=tc[2];
               ec3:(A.c)[3][1]=tc[3];

               afsol:solve([ea1,ea2,ea3,eb1,eb2,eb3,ec1,ec2,ec3],
               [a11,a12,a13,a21,a22,a23,a31,a32,a33]);

               gives: (with a11=A(0,0), ... )
            */
            double a,b,c,d,
              ad,bc,dlt,ap,bp,cp,dp,aa,bb,ab,cc,cd,dd;
            a=xb-xa;
            b=xc-xa;
            c=yb-ya;
            d=yc-ya;
            ad=a*d;
            bc=b*c;
            dlt=ad-bc;
            ap=d/dlt;
            bp=-b/dlt;
            cp=-c/dlt;
            dp=a/dlt;
            aa=a*a;
            bb=b*b;
            ab=a*b;
            cc=c*c;
            cd=c*d;
            dd=d*d;


	    
	    A(0,0) =  ap;
	    A(0,1) =  bp;
	    A(0,2) =  -ap*xa-bp*ya;
	    A(1,0) =  cp;
	    A(1,1) =  dp;
	    A(1,2) =  -cp*xa-dp*ya;
	    A(2,0) =  0.0;
	    A(2,1) =  0.0;
	    A(2,2) =  1.0;


            for(int l=0;l<3;l++)
              for(int m=0;m<3;m++)
                AfTr(3*i+o,3*m+l) = A(l,m);

            /*
            // test:
            */
            VectorXd ta(3),tb(3),tc(3);
            VectorXd xya(3),xyb(3),xyc(3);
            xya << xa, ya, 1.0;
            xyb << xb, yb, 1.0;
            xyc << xc, yc, 1.0;
            ta=A*xya;
            tb=A*xyb;
            tc=A*xyc;
            /*
              Test passed, ok.
            */

            //
            double z_a,z_b,z_c,
              z_x_a,z_x_b,z_x_c,
              z_y_a,z_y_b,z_y_c,
              z_xy_a,z_xy_b,z_xy_c,
              z_xx_a,z_xx_b,z_xx_c,
              z_yy_a,z_yy_b,z_yy_c;


            z_a=zT[0];
            z_b=zT[1];
            z_c=zT[2];
            z_x_a=Z_x[iT[0]];
            z_x_b=Z_x[iT[1]];
            z_x_c=Z_x[iT[2]];
            z_y_a=Z_y[iT[0]];
            z_y_b=Z_y[iT[1]];
            z_y_c=Z_y[iT[2]];
            z_xy_a=Z_xy[iT[0]];
            z_xy_b=Z_xy[iT[1]];
            z_xy_c=Z_xy[iT[2]];
            z_xx_a=Z_xx[iT[0]];
            z_xx_b=Z_xx[iT[1]];
            z_xx_c=Z_xx[iT[2]];
            z_yy_a=Z_yy[iT[0]];
            z_yy_b=Z_yy[iT[1]];
            z_yy_c=Z_yy[iT[2]];


            double z_u_a,z_u_b,z_u_c,
              z_v_a,z_v_b,z_v_c,
              z_uv_a,z_uv_b,z_uv_c,
              z_uu_a,z_uu_b,z_uu_c,
              z_vv_a,z_vv_b,z_vv_c;
	      // z_uv_ab, z_uv_bc,z_uv_ca;


            // gradients and Hesse matrices in x-y and u-v coordinates:
            VectorXd grad_xy_a(2), grad_xy_b(2), grad_xy_c(2);
            VectorXd grad_uv_a(2), grad_uv_b(2), grad_uv_c(2);
            MatrixXd H_xy_a(2,2), H_uv_a(2,2),
              H_xy_b(2,2), H_uv_b(2,2),
              H_xy_c(2,2), H_uv_c(2,2);

            grad_xy_a << z_x_a, z_y_a;
            grad_xy_b << z_x_b, z_y_b;
            grad_xy_c << z_x_c, z_y_c;

            // transform gradients, take only the stretch/rotation part of A

            grad_uv_a = A.block(0,0,2,2).inverse().transpose()*grad_xy_a;
            grad_uv_b = A.block(0,0,2,2).inverse().transpose()*grad_xy_b;
            grad_uv_c = A.block(0,0,2,2).inverse().transpose()*grad_xy_c;


            H_xy_a << z_xx_a, z_xy_a,
              z_xy_a, z_yy_a;
            H_xy_b << z_xx_b, z_xy_b,
              z_xy_b, z_yy_b;
            H_xy_c << z_xx_c, z_xy_c,
              z_xy_c, z_yy_c;

            // transform Hesse matrices:
            H_uv_a = A.block(0,0,2,2).inverse()*H_xy_a*A.block(0,0,2,2).inverse().transpose();
            H_uv_b = A.block(0,0,2,2).inverse()*H_xy_b*A.block(0,0,2,2).inverse().transpose();
            H_uv_c = A.block(0,0,2,2).inverse()*H_xy_c*A.block(0,0,2,2).inverse().transpose();

            // directional derivatives:
            // vectors forming the triangle :
            VectorXd ab_xy(2),bc_xy(2),ca_xy(2),n_xy_ab(2),n_xy_bc(2),n_xy_ca(2);
            ab_xy << xb-xa, yb-ya;
            bc_xy << xc-xb, yc-yb;
            ca_xy << xa-xc, ya-yc;
            MatrixXd N(2,2);
            N << 0.0, 1.0, -1.0, 0.0;


            // make normal vectors:
            VectorXd n_uv_ab(2), n_uv_bc(2), n_uv_ca(2);
            n_xy_ab = N*ab_xy;
            n_xy_ab = n_xy_ab / sqrt(n_xy_ab.norm());
            n_xy_bc = N*bc_xy;
            n_xy_bc = n_xy_bc / sqrt(n_xy_bc.norm());
            n_xy_ca = N*ca_xy;
            n_xy_ca = n_xy_ca / sqrt(n_xy_ca.norm());
            // after affine transformation
            VectorXd ab_uv(2),bc_uv(2),ca_uv(2);
            //ab_uv << 0.0, -1.0;
            ab_uv << 1.0, 0.0;
            bc_uv << -1.0, 1.0;
            //ca_uv << 1.0, 0.0;
            ca_uv << 0.0, -1.0;
            n_uv_ab = N*ab_uv;
            n_uv_ab = n_uv_ab / sqrt(n_uv_ab.norm());
            n_uv_bc = N*bc_uv;
            n_uv_bc = n_uv_bc / sqrt(n_uv_bc.norm());
            n_uv_ca = N*ca_uv;
            n_uv_ca = n_uv_ca / sqrt(n_uv_ca.norm());

            // gradients at midpoints
            VectorXd grad_xy_ab(2), grad_xy_bc(2), grad_xy_ca(2);
            grad_xy_ab << Z_x_ab[i], Z_y_ab[i];
            grad_xy_bc << Z_x_bc[i], Z_y_bc[i];
            grad_xy_ca << Z_x_ca[i], Z_y_ca[i];

            // transform to uv
            VectorXd grad_uv_ab(2), grad_uv_bc(2), grad_uv_ca(2);
            grad_uv_ab=A.block(0,0,2,2).inverse().transpose()*grad_xy_ab;
            grad_uv_bc=A.block(0,0,2,2).inverse().transpose()*grad_xy_bc;
            grad_uv_ca=A.block(0,0,2,2).inverse().transpose()*grad_xy_ca;

            // directional derivative in uv coords:
            // FIXME: + or - ?
            // z_uv_ab = grad_uv_ab.transpose()*n_uv_ab;
            // z_uv_bc = grad_uv_bc.transpose()*n_uv_bc;
            // z_uv_ca = grad_uv_ca.transpose()*n_uv_ca;

            // extract transformed values:
            z_u_a=grad_uv_a[0];
            z_v_a=grad_uv_a[1];
            z_u_b=grad_uv_b[0];
            z_v_b=grad_uv_b[1];
            z_u_c=grad_uv_c[0];
            z_v_c=grad_uv_c[1];

            z_uu_a=H_uv_a(0,0);
            z_uv_a=H_uv_a(1,0);
            z_vv_a=H_uv_a(1,1);

            z_uu_b=H_uv_b(0,0);
            z_uv_b=H_uv_b(1,0);
            z_vv_b=H_uv_b(1,1);

            z_uu_c=H_uv_c(0,0);
            z_uv_c=H_uv_c(1,0);
            z_vv_c=H_uv_c(1,1);



	    a00 = z_a;
	    a10 = z_u_a;
	    a01 = z_v_a;
	    a20 = z_uu_a*0.5;
	    a11 = z_uv_a;
	    a02 = z_vv_a*0.5;

	    double h1 = z_b-a00-a10-a20;
	    double h2 = z_u_b-a10-z_uu_a;
	    double h3 = z_uu_b-z_uu_a;
	    a30 = 10.0*h1-4.0*h2+0.5*h3;
	    a40 = -15.0*h1+7.0*h2-h3;
	    a50 = 6.0*h1-3.0*h2+0.5*h3;

	    h1 = z_c - a00 - a01 - a02;
	    h2 = z_v_c - a01 - z_vv_a;
	    h3 = z_vv_c - z_vv_a;
	    a03 = 10.0*h1-4.0*h2+0.5*h3;
	    a04 = -15.0*h1+7.0*h2-h3;
	    a05 = 6.0*h1-3.0*h2+0.5*h3;

	    double lusq=aa+cc;
	    double lvsq=bb+dd;
	    double spuv=ab+cd;
	    a41 = 5.0*spuv/lusq*a50;
	    a14 = 5.0*spuv/lvsq*a05;

	    h1 = z_v_b-a01-a11-a41;
	    h2 = z_uv_b -a11-4.0*a41;
	    a21 = 3.0*h1-h2;
	    a31 = -2.0*h1+h2;

	    h1 = z_u_c-a10-a11-a14;
	    h2 = z_uv_c -a11-4.0*a14;
	    a12 = 3.0*h1-h2;
	    a13 = -2.0*h1+h2;

	    double e1=(lvsq-spuv)/((lvsq-spuv)+(lusq-spuv));
	    double e2=1.0-e1;
	    double g1=5.0*e1-2.0;
	    double g2=1.0-g1;
	    h1 = 5.0*(e1*(a50-a41)+e2*(a05-a14))+(a41+a14);
	    h2 = 0.5* z_vv_b -a02-a12;
	    h3 = 0.5* z_uu_c -a20-a21;
	    a22 = h1+g1*h2+g2*h3;
	    a32 = h2-a22;
	    a23 = h3-a22;



            // TODO: compare estimates to polynomial derivatives
            // if differences are too large fall back to lower degree in pD
            // per triangle, use transformed coordinates

            // double p_a, p_b, p_c,     // polynom values
            //   p_u_a, p_u_b, p_u_c,    // first derivatives
            //   p_v_a, p_v_b, p_v_c,
            //   p_uu_a, p_uv_a, p_vv_a, // second -"-
            //   p_uu_b, p_uv_b, p_vv_b,
            //   p_uu_c, p_uv_c, p_vv_c,
            //   p_ab, p_bc, p_ca; // directional derivatives in midpoints
            //      in u,v coordinates and counterclockwise it should hold:
            //                    p_ab=-p_v_a,
            //                    p_bc=sqrt(2)/2(p_u_a+p_v_a)
            //                    p_ca=-p_u_a,
            // double ujk,vjk;
            // ujk=0.0;vjk=0.0;


            /*
              (%i77) factorsum(subst([x=ujk,y=vjk],f(x,y)));
              (%o77) vjk (a01 + ujk (ujk (a21 + ujk (a41 ujk + a31)) + a11)
              + vjk (vjk (a03 + ujk (a23 ujk + a13) + vjk (a05 vjk + a14 ujk + a04))
              + ujk (a12 + ujk (a32 ujk + a22)) + a02))
              + ujk (a10 + ujk (ujk (a30 + ujk (a50 ujk + a40)) + a20)) + a00
            */
            // p_a=vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00;
            /*
              (%i74) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,1)));
              (%o74) vjk (a11 + ujk (ujk (3 a31 + 4 a41 ujk) + 2 a21)
              + vjk (vjk (a13 + 2 a23 ujk + a14 vjk) + ujk (2 a22 + 3 a32 ujk) + a12))
              + ujk (2 a20 + ujk (ujk (4 a40 + 5 a50 ujk) + 3 a30)) + a10
            */

            // p_u_a=vjk*(a11+ujk*(ujk*(3.0*a31+4*a41*ujk)+2.0*a21)+vjk*(vjk*(a13+2.0*a23*ujk+a14*vjk)+ujk*(2.0*a22+3.0*a32*ujk)+a12))+ujk*(2.0*a20+ujk*(ujk*(4.0*a40+5.0*a50*ujk)+3.0*a30))+a10;

            // p_v_a=vjk*(2.0*a02+2*ujk*(ujk*(a22+a32*ujk)+a12)+vjk*(vjk*(4.0*(a14*ujk+a04)+5.0*a05*vjk)+3.0*ujk*(a13+a23*ujk)+3.0*a03))+ujk*(a11+ujk*(ujk*(a31+a41*ujk)+a21))+a01;


            /*
              (%i84) factorsum(subst([x=ujk,y=vjk],diff(diff(f(x,y),x,1),y,1)));
              (%o84) vjk (2 a12 + 2 ujk (3 a32 ujk + 2 a22)
              + vjk (4 a14 vjk + 6 a23 ujk + 3 a13)) + ujk (2 a21 + ujk (4 a41 ujk + 3 a31))
              + a11
            */

            // p_uv_a=vjk*(2.0*a12+2.0*ujk*(3.0*a32*ujk+2.0*a22)+vjk*(4.0*a14*vjk+6.0*a23*ujk+3.0*a13))+ujk*(2.0*a21+ujk*(4.0*a41*ujk+3.0*a31))+a11;

            /*
              (%i80) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,2)));
              (%o80) 2 (vjk (a21 + 3 ujk (2 a41 ujk + a31) + vjk (a23 vjk + 3 a32 ujk + a22))
              + ujk (3 a30 + 2 ujk (5 a50 ujk + 3 a40)) + a20)

            */
            // p_uu_a=2.0*(vjk*(a21+3.0*ujk*(2.0*a41*ujk+a31)+vjk*(a23*vjk+3.0*a32*ujk+a22))+ujk*(3.0*a30+2*ujk*(5.0*a50*ujk+3.0*a40))+a20);

            /*
              (%i81) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),y,2)));
              (%o81) 2 (vjk (3 a03 + 3 ujk (a23 ujk + a13)
              + 2 vjk (5 a05 vjk + 3 (a04 + a14 ujk))) + ujk (a12 + ujk (a32 ujk + a22))
              + a02)
            */
            // p_vv_a=2.0*(vjk*(3.0*a03+3.0*ujk*(a23*ujk+a13)+2.0*vjk*(5.0*a05*vjk+3.0*(a04+a14*ujk)))+ujk*(a12+ujk*(a32*ujk+a22))+a02);

            /*
              string(factorsum(subst([x=0.5,y=0],diff(f(x,y),y,1))));

              (%o6)                  (a41+2*a31+4*a21+8*a11+16*a01)/16
            */
            //p_ab=-(a41+2.0*a31+4.0*a21+8.0*a11+16.0*a01)/16.0;


            // ujk=0.0;vjk=1.0;
            /*
              (%i77) factorsum(subst([x=ujk,y=vjk],f(x,y)));
              (%o77) vjk (a01 + ujk (ujk (a21 + ujk (a41 ujk + a31)) + a11)
              + vjk (vjk (a03 + ujk (a23 ujk + a13) + vjk (a05 vjk + a14 ujk + a04))
              + ujk (a12 + ujk (a32 ujk + a22)) + a02))
              + ujk (a10 + ujk (ujk (a30 + ujk (a50 ujk + a40)) + a20)) + a00
            */
            // p_b=vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00;
            /*
              (%i74) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,1)));
              (%o74) vjk (a11 + ujk (ujk (3 a31 + 4 a41 ujk) + 2 a21)
              + vjk (vjk (a13 + 2 a23 ujk + a14 vjk) + ujk (2 a22 + 3 a32 ujk) + a12))
              + ujk (2 a20 + ujk (ujk (4 a40 + 5 a50 ujk) + 3 a30)) + a10
            */

            // p_u_b=vjk*(a11+ujk*(ujk*(3.0*a31+4*a41*ujk)+2.0*a21)+vjk*(vjk*(a13+2.0*a23*ujk+a14*vjk)+ujk*(2.0*a22+3.0*a32*ujk)+a12))+ujk*(2.0*a20+ujk*(ujk*(4.0*a40+5.0*a50*ujk)+3.0*a30))+a10;

            // p_v_b=vjk*(2.0*a02+2*ujk*(ujk*(a22+a32*ujk)+a12)+vjk*(vjk*(4.0*(a14*ujk+a04)+5.0*a05*vjk)+3.0*ujk*(a13+a23*ujk)+3.0*a03))+ujk*(a11+ujk*(ujk*(a31+a41*ujk)+a21))+a01;


            /*
              (%i84) factorsum(subst([x=ujk,y=vjk],diff(diff(f(x,y),x,1),y,1)));
              (%o84) vjk (2 a12 + 2 ujk (3 a32 ujk + 2 a22)
              + vjk (4 a14 vjk + 6 a23 ujk + 3 a13)) + ujk (2 a21 + ujk (4 a41 ujk + 3 a31))
              + a11
            */

            // p_uv_b=vjk*(2.0*a12+2.0*ujk*(3.0*a32*ujk+2.0*a22)+vjk*(4.0*a14*vjk+6.0*a23*ujk+3.0*a13))+ujk*(2.0*a21+ujk*(4.0*a41*ujk+3.0*a31))+a11;

            /*
              (%i80) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,2)));
              (%o80) 2 (vjk (a21 + 3 ujk (2 a41 ujk + a31) + vjk (a23 vjk + 3 a32 ujk + a22))
              + ujk (3 a30 + 2 ujk (5 a50 ujk + 3 a40)) + a20)

            */
            // p_uu_b=2.0*(vjk*(a21+3.0*ujk*(2.0*a41*ujk+a31)+vjk*(a23*vjk+3.0*a32*ujk+a22))+ujk*(3.0*a30+2*ujk*(5.0*a50*ujk+3.0*a40))+a20);

            /*
              (%i81) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),y,2)));
              (%o81) 2 (vjk (3 a03 + 3 ujk (a23 ujk + a13)
              + 2 vjk (5 a05 vjk + 3 (a04 + a14 ujk))) + ujk (a12 + ujk (a32 ujk + a22))
              + a02)
            */
            // p_vv_b=2.0*(vjk*(3.0*a03+3.0*ujk*(a23*ujk+a13)+2.0*vjk*(5.0*a05*vjk+3.0*(a04+a14*ujk)))+ujk*(a12+ujk*(a32*ujk+a22))+a02);

            /*
              factorsum(subst([x=0.5,y=0.5],sqrt(2)/2*(diff(f(x,y),x,1)+diff(f(x,y),y,1))));
              (5*(a05+a14+a23+a32+a41+a50)+8*(a04+a13+a22+a31+a40)+12*(a03+a12+a21+a30 \
              )+16*(a01+a02+a10+a11+a20))/2^(9/2)
            */
            //double p_2_9_2=pow(2.0,4.5);
            //p_bc=(5.0*(a05+a14+a23+a32+a41+a50)+8.0*(a04+a13+a22+a31+a40)+12.0*(a03+a12+a21+a30)+16.0*(a01+a02+a10+a11+a20))/p_2_9_2;

            // ujk=1.0;vjk=0.0;

            /*
              (%i77) factorsum(subst([x=ujk,y=vjk],f(x,y)));
              (%o77) vjk (a01 + ujk (ujk (a21 + ujk (a41 ujk + a31)) + a11)
              + vjk (vjk (a03 + ujk (a23 ujk + a13) + vjk (a05 vjk + a14 ujk + a04))
              + ujk (a12 + ujk (a32 ujk + a22)) + a02))
              + ujk (a10 + ujk (ujk (a30 + ujk (a50 ujk + a40)) + a20)) + a00
            */
            // p_c=vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00;
            /*
              (%i74) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,1)));
              (%o74) vjk (a11 + ujk (ujk (3 a31 + 4 a41 ujk) + 2 a21)
              + vjk (vjk (a13 + 2 a23 ujk + a14 vjk) + ujk (2 a22 + 3 a32 ujk) + a12))
              + ujk (2 a20 + ujk (ujk (4 a40 + 5 a50 ujk) + 3 a30)) + a10
            */

            // p_u_c=vjk*(a11+ujk*(ujk*(3.0*a31+4*a41*ujk)+2.0*a21)+vjk*(vjk*(a13+2.0*a23*ujk+a14*vjk)+ujk*(2.0*a22+3.0*a32*ujk)+a12))+ujk*(2.0*a20+ujk*(ujk*(4.0*a40+5.0*a50*ujk)+3.0*a30))+a10;

            // p_v_c=vjk*(2.0*a02+2*ujk*(ujk*(a22+a32*ujk)+a12)+vjk*(vjk*(4.0*(a14*ujk+a04)+5.0*a05*vjk)+3.0*ujk*(a13+a23*ujk)+3.0*a03))+ujk*(a11+ujk*(ujk*(a31+a41*ujk)+a21))+a01;


            /*
              (%i84) factorsum(subst([x=ujk,y=vjk],diff(diff(f(x,y),x,1),y,1)));
              (%o84) vjk (2 a12 + 2 ujk (3 a32 ujk + 2 a22)
              + vjk (4 a14 vjk + 6 a23 ujk + 3 a13)) + ujk (2 a21 + ujk (4 a41 ujk + 3 a31))
              + a11
            */

            // p_uv_c=vjk*(2.0*a12+2.0*ujk*(3.0*a32*ujk+2.0*a22)+vjk*(4.0*a14*vjk+6.0*a23*ujk+3.0*a13))+ujk*(2.0*a21+ujk*(4.0*a41*ujk+3.0*a31))+a11;

            /*
              (%i80) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),x,2)));
              (%o80) 2 (vjk (a21 + 3 ujk (2 a41 ujk + a31) + vjk (a23 vjk + 3 a32 ujk + a22))
              + ujk (3 a30 + 2 ujk (5 a50 ujk + 3 a40)) + a20)

            */
            // p_uu_c=2.0*(vjk*(a21+3.0*ujk*(2.0*a41*ujk+a31)+vjk*(a23*vjk+3.0*a32*ujk+a22))+ujk*(3.0*a30+2*ujk*(5.0*a50*ujk+3.0*a40))+a20);

            /*
              (%i81) factorsum(subst([x=ujk,y=vjk],diff(f(x,y),y,2)));
              (%o81) 2 (vjk (3 a03 + 3 ujk (a23 ujk + a13)
              + 2 vjk (5 a05 vjk + 3 (a04 + a14 ujk))) + ujk (a12 + ujk (a32 ujk + a22))
              + a02)
            */
            // p_vv_c=2.0*(vjk*(3.0*a03+3.0*ujk*(a23*ujk+a13)+2.0*vjk*(5.0*a05*vjk+3.0*(a04+a14*ujk)))+ujk*(a12+ujk*(a32*ujk+a22))+a02);

            /*
              string(factorsum(subst([x=0,y=0.5],diff(f(x,y),x,1))));

              (%o9)                  (a14+2*a13+4*a12+8*a11+16*a10)/16
            */

            //p_ca=-(a14+2.0*a13+4.0*a12+8.0*a11+16.0*a10)/16.0;

            /*
            double alpha_1=1.0; // ? downweight errors in first deriv.
            double alpha_2=1.0; // ? downweight errors in second ..
            double err_a=(z_a-p_a)*(z_a-p_a)+
              alpha_1*((z_u_a-p_u_a)*(z_u_a-p_u_a)+(z_v_a-p_v_a)*(z_v_a-p_v_a))+
              alpha_2*((z_uu_a-p_uu_a)*(z_uu_a-p_uu_a)+(z_uv_a-p_uv_a)*(z_uv_a-p_uv_a)+(z_vv_a-p_vv_a)*(z_vv_a-p_vv_a));
            double err_b=(z_b-p_b)*(z_b-p_b)+
              alpha_1*((z_u_b-p_u_b)*(z_u_b-p_u_b)+(z_v_b-p_v_b)*(z_v_b-p_v_b))+
              alpha_2*((z_uu_b-p_uu_b)*(z_uu_b-p_uu_b)+(z_uv_b-p_uv_b)*(z_uv_b-p_uv_b)+(z_vv_b-p_vv_b)*(z_vv_b-p_vv_b));
            double err_c=(z_c-p_c)*(z_c-p_c)+
              alpha_1*((z_u_c-p_u_c)*(z_u_c-p_u_c)+(z_v_c-p_v_c)*(z_v_c-p_v_c))+
              alpha_2*((z_uu_c-p_uu_c)*(z_uu_c-p_uu_c)+(z_uv_c-p_uv_c)*(z_uv_c-p_uv_c)+(z_vv_c-p_vv_c)*(z_vv_c-p_vv_c));
            */
	    // FIXME: continue to use err_a, err_b, err_c 
	    
            par_o << a00,a01,a02,a03,a04,a05,a10,a11,a12,a13,a14,a20,a21,a22,a23,a30,a31,a32,a40,a41,a50;
            for(int k=0;k<21;k++){
              par(3*i+o,k) =par_o[k];
            }
          } // end loop over o
        } // iterate over triangles (i)

        for(int e=0;e<EdgeCheck.nE;e++){
          for(int l=0;l<nCheck; l++){
            EdgeCheck.zBl(e,l)=0;//std::nan("1");
            EdgeCheck.zBr(e,l)=0;//std::nan("1");
          }
        }

        for(int e=0;e<EdgeCheck.nE;e++){
          //Rcout << "edge " << e+1 << " is in triangle " << EdgeCheck.t1[e] << " and " << EdgeCheck.t2[e] << std::endl ;
          //Rcout << "with nodes (" << xD[EdgeCheck.i1[e]] << ", " << yD[EdgeCheck.i1[e]] <<
          //  ") and ("  << xD[EdgeCheck.i2[e]] << ", " << yD[EdgeCheck.i2[e]] << ")" << std::endl;
          //Rcout << "test points are (" << EdgeCheck.xB(e,0) << ", " <<  EdgeCheck.yB(e,0) << ") (" << EdgeCheck.xB(e,1) << ", " <<  EdgeCheck.yB(e,1) << ") (" << EdgeCheck.xB(e,2) << ", " <<  EdgeCheck.yB(e,2) << ") " << std::endl;

          for(int ik=0;ik<2;ik++){
            int k;
            if(ik==0)
              k=EdgeCheck.t1[e];
            else
              k=EdgeCheck.t2[e];
            if(k>0){
              //Rcout << "e: " << e << " k: " << k << std::endl;
              int barycycles=1;
              if(baryweight)
                barycycles=3;
              for(int o=0;o<barycycles;o++){
                if(o==0){
                  iT[0]=tXY.i1[k]; iT[1]=tXY.i2[k]; iT[2]=tXY.i3[k];
                  xT[0]=xD[tXY.i1[k]]; xT[1]=xD[tXY.i2[k]]; xT[2]=xD[tXY.i3[k]];
                  yT[0]=yD[tXY.i1[k]]; yT[1]=yD[tXY.i2[k]]; yT[2]=yD[tXY.i3[k]];
                  zT[0]=zD[tXY.i1[k]]; zT[1]=zD[tXY.i2[k]]; zT[2]=zD[tXY.i3[k]];
                } else if(o==1){
                  iT[0]=tXY.i2[k]; iT[1]=tXY.i3[k]; iT[2]=tXY.i1[k];
                  xT[0]=xD[tXY.i2[k]]; xT[1]=xD[tXY.i3[k]]; xT[2]=xD[tXY.i1[k]];
                  yT[0]=yD[tXY.i2[k]]; yT[1]=yD[tXY.i3[k]]; yT[2]=yD[tXY.i1[k]];
                  zT[0]=zD[tXY.i2[k]]; zT[1]=zD[tXY.i3[k]]; zT[2]=zD[tXY.i1[k]];
                } else {
                  iT[0]=tXY.i3[k]; iT[1]=tXY.i1[k]; iT[2]=tXY.i2[k];
                  xT[0]=xD[tXY.i3[k]]; xT[1]=xD[tXY.i1[k]]; xT[2]=xD[tXY.i2[k]];
                  yT[0]=yD[tXY.i3[k]]; yT[1]=yD[tXY.i1[k]]; yT[2]=yD[tXY.i2[k]];
                  zT[0]=zD[tXY.i3[k]]; zT[1]=zD[tXY.i1[k]]; zT[2]=zD[tXY.i2[k]];
                }

                a00=par(3*k+o,0);
                a01=par(3*k+o,1);
                a02=par(3*k+o,2);
                a03=par(3*k+o,3);
                a04=par(3*k+o,4);
                a05=par(3*k+o,5);
                a10=par(3*k+o,6);
                a11=par(3*k+o,7);
                a12=par(3*k+o,8);
                a13=par(3*k+o,9);
                a14=par(3*k+o,10);
                a20=par(3*k+o,11);
                a21=par(3*k+o,12);
                a22=par(3*k+o,13);
                a23=par(3*k+o,14);
                a30=par(3*k+o,15);
                a31=par(3*k+o,16);
                a32=par(3*k+o,17);
                a40=par(3*k+o,18);
                a41=par(3*k+o,19);
                a50=par(3*k+o,20);

                A(0,0) =  AfTr(3*k+o,0);
                A(1,0) =  AfTr(3*k+o,1);
                A(2,0) =  AfTr(3*k+o,2);
                A(0,1) =  AfTr(3*k+o,3);
                A(1,1) =  AfTr(3*k+o,4);
                A(2,1) =  AfTr(3*k+o,5);
                A(0,2) =  AfTr(3*k+o,6);
                A(1,2) =  AfTr(3*k+o,7);
                A(2,2) =  AfTr(3*k+o,8);
                for(int l=0;l<3;l++){
                  // barycentric:
                  double a = ((yT[1] - yT[2])*(EdgeCheck.xB(e,l) - xT[2]) + (xT[2] - xT[1])*(EdgeCheck.yB(e,l) - yT[2])) /
                    ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));

                  //double b = ((yT[2] - yT[0])*(EdgeCheck.xB(e,l) - xT[2]) + (xT[0] - xT[2])*(EdgeCheck.yB(e,l) - yT[2])) /
                  //  ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));

                  // double c = 1 - a - b;
                  VectorXd xy(3), uv(3);

                  // use homogeneous coordinates:
                  xy << EdgeCheck.xB(e,l), EdgeCheck.yB(e,l), 1.0;
                  double ujk,vjk;

                  uv = A*xy;
                  ujk=uv[0];vjk=uv[1];

                  if(!baryweight)
                    a=1.0;

                  if(o==0){
                    if(ik==0)
                      EdgeCheck.zBl(e,l)= a*(vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00 ) ;
                    else
                      EdgeCheck.zBr(e,l)= a*(vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00 ) ;
                  } else {
                    if(ik==0)
                      EdgeCheck.zBl(e,l)+= a*( vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00  );
                    else
                      EdgeCheck.zBr(e,l)+= a*( vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00  );
                  }
                } // l
              } // o
            } // k>0
          } // zBl, zBr

        }
        //Rcout << "Edge check: left " << EdgeCheck.zBl << std::endl;
        //Rcout << "Edge check: right " << EdgeCheck.zBr << std::endl;
        VectorXd EdgeDelta=((EdgeCheck.zBl-EdgeCheck.zBr).array()*(EdgeCheck.zBl-EdgeCheck.zBr).array()).rowwise().sum();

        if(autodegree){
          for(int e=0;e<EdgeCheck.nE;e++){
            //Rcout << "sum delta^2(" << e << ")=" << EdgeDelta[e] << std::endl;
            if(!std::isinf(EdgeDelta[e])){
              if(EdgeDelta[e]>adtol){
                //Rcout << "reducing degree for node " << EdgeCheck.i1[e] << " and " << EdgeCheck.i2[e] << std::endl;

                // dont fall back to degree=0 !!!
                if(doEstD[EdgeCheck.i1[e]]>1){
                  doEstD[EdgeCheck.i1[e]]=doEstD[EdgeCheck.i1[e]]-1 ;
                  //Rcout << " to " <<  doEstD[EdgeCheck.i1[e]] << std::endl;
                }
                if(doEstD[EdgeCheck.i2[e]]>1){
                  doEstD[EdgeCheck.i2[e]]=doEstD[EdgeCheck.i2[e]]-1 ;
                  //Rcout << " to " <<  doEstD[EdgeCheck.i2[e]] << std::endl;
                }
              }
            }
          }
        }
      } // end for dg
      //for(int i=0; i<nxD; i++)
      //  Rcout << "data point " << i << " is estimated with degree " << doEstD[i] << std::endl;
    } // !linear
   
    //###################################### WORK END
    
    
    
    // part 2
    // iterate over triangles
    for(int i=0;i<nT;i++){
      //################################### WORK 2 START
      if(!linear){
        int barycycles=1;
        if(baryweight)
          barycycles=3;
        for(int i_barycycle=0;i_barycycle<barycycles;i_barycycle++){



          if(i_barycycle==0){
            iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
            xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
            yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];
            zT[0]=zD[tXY.i1[i]]; zT[1]=zD[tXY.i2[i]]; zT[2]=zD[tXY.i3[i]];
          } else if(i_barycycle==1){
            iT[0]=tXY.i2[i]; iT[1]=tXY.i3[i]; iT[2]=tXY.i1[i];
            xT[0]=xD[tXY.i2[i]]; xT[1]=xD[tXY.i3[i]]; xT[2]=xD[tXY.i1[i]];
            yT[0]=yD[tXY.i2[i]]; yT[1]=yD[tXY.i3[i]]; yT[2]=yD[tXY.i1[i]];
            zT[0]=zD[tXY.i2[i]]; zT[1]=zD[tXY.i3[i]]; zT[2]=zD[tXY.i1[i]];
          } else {
            iT[0]=tXY.i3[i]; iT[1]=tXY.i1[i]; iT[2]=tXY.i2[i];
            xT[0]=xD[tXY.i3[i]]; xT[1]=xD[tXY.i1[i]]; xT[2]=xD[tXY.i2[i]];
            yT[0]=yD[tXY.i3[i]]; yT[1]=yD[tXY.i1[i]]; yT[2]=yD[tXY.i2[i]];
            zT[0]=zD[tXY.i3[i]]; zT[1]=zD[tXY.i1[i]]; zT[2]=zD[tXY.i2[i]];
          }
          /*
            if(orientation(i,i_barycycle)==0.0)
            Rf_error("triangle collapsed!");
            // prefer counter clockwise orientation:
            if(orientation(i,i_barycycle)>0.0){
            // swap points 1 and 2
            double swap_x, swap_y, swap_z;
            int swap_i;

            swap_i=iT[2]; iT[2]=iT[1]; iT[1]=swap_i;
            swap_x=xT[2]; xT[2]=xT[1]; xT[1]=swap_x;
            swap_y=yT[2]; yT[2]=yT[1]; yT[1]=swap_y;
            swap_z=zT[2]; zT[2]=zT[1]; zT[1]=swap_z;

            }
          */

          a00=par(3*i+i_barycycle,0);
          a01=par(3*i+i_barycycle,1);
          a02=par(3*i+i_barycycle,2);
          a03=par(3*i+i_barycycle,3);
          a04=par(3*i+i_barycycle,4);
          a05=par(3*i+i_barycycle,5);
          a10=par(3*i+i_barycycle,6);
          a11=par(3*i+i_barycycle,7);
          a12=par(3*i+i_barycycle,8);
          a13=par(3*i+i_barycycle,9);
          a14=par(3*i+i_barycycle,10);
          a20=par(3*i+i_barycycle,11);
          a21=par(3*i+i_barycycle,12);
          a22=par(3*i+i_barycycle,13);
          a23=par(3*i+i_barycycle,14);
          a30=par(3*i+i_barycycle,15);
          a31=par(3*i+i_barycycle,16);
          a32=par(3*i+i_barycycle,17);
          a40=par(3*i+i_barycycle,18);
          a41=par(3*i+i_barycycle,19);
          a50=par(3*i+i_barycycle,20);


          A(0,0) =  AfTr(3*i+i_barycycle,0);
          A(1,0) =  AfTr(3*i+i_barycycle,1);
          A(2,0) =  AfTr(3*i+i_barycycle,2);
          A(0,1) =  AfTr(3*i+i_barycycle,3);
          A(1,1) =  AfTr(3*i+i_barycycle,4);
          A(2,1) =  AfTr(3*i+i_barycycle,5);
          A(0,2) =  AfTr(3*i+i_barycycle,6);
          A(1,2) =  AfTr(3*i+i_barycycle,7);
          A(2,2) =  AfTr(3*i+i_barycycle,8);



          if(as<std::string>(output)=="grid"){
            // iterate only over grid points (j,k) inside bounding box of triangle i
            for(int j=jTsw[i]; j<jTne[i]; j++) {
              for(int k=kTsw[i]; k<kTne[i]; k++) {
                // calculate barycentric coordinates:
                double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[k] - yT[2])) /
                  ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
                double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[k] - yT[2])) /
                  ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
                double c = 1 - a - b;
                // check if inside triangle, handle only yet untouched grid points
                //if(R_IsNA(z(j,k))){
                if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
                  if(linear)
                    z(j,k)=a*zT[0]+b*zT[1]+c*zT[2]; // FIXME: reached? duplicate below!!
                  else{
                    // affine transformation
                    VectorXd xy(3), uv(3);
                    // use homogeneous coordinates:
                    xy << x[j], y[k], 1.0;
                    double ujk,vjk;
                    
                    uv = A*xy;
                    ujk=uv[0];vjk=uv[1];
                    
                    /*
                      (%i72) factorsum(f(x,y));
                      (%o72) y (a01 + x (x (a21 + x (a41 x + a31)) + a11)
                      + y (y (a03 + x (a23 x + a13) + y (a05 y + a14 x + a04))
                      + x (a12 + x (a32 x + a22)) + a02)) + x
                      (a10 + x (x (a30 + x (a50 x + a40)) + a20)) + a00
                    */
                    // Use barycentric coordinates as weights for
                    // convex linear combination of the three polynomials:
                    // variable i_barycycle iterates over the triangle enumerations, so if a corner
                    // is named "a" for i_barycycle==0 it will be "b" for i_barycycle==1 and "c" for i_barycycle==3
                    // So the sum of all barycentric "a"s over i_barycycle=0,1,2 will give 1
                    // as they change their meaning relative to the first enumeration of the triangle.
                    // Corner "a" will be transformed to (u,v)=(0,0) and will have the smallest errors.
                    // For that reason it seems ok to use the barycentric coordniates (a,b,c) of the pixel
                    // (ujk,vjik) as convex linear combination weights. A pixel near one of the corners "a"
                    // "b" or "c" will then have minimal error. The barycentric coordinates (a,b,c) can
                    // also be expressed as (a,a',a'') where (ujk,vjik)=(a,b,c) in enumeration i_barycycle==0,
                    // (ujk,vjik)=(a',b',c') in enumeration i_barycycle==1 and (ujk,vjik)=(a'',b'',c'') for i_barycycle==2
                    
                    // no baryweight: only single term in summation
                    if(!baryweight)
                      a=1.0;
                    
                    if(i_barycycle==0){
                      z(j,k)=a*(  vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00  );
                    } else {
                      z(j,k)+=a*(  vjk*(a01+ujk*(ujk*(a21+ujk*(a41*ujk+a31))+a11)+vjk*(vjk*(a03+ujk*(a23*ujk+a13)+vjk*(a05*vjk+a14*ujk+a04))+ujk*(a12+ujk*(a32*ujk+a22))+a02))+ujk*(a10+ujk*(ujk*(a30+ujk*(a50*ujk+a40))+a20))+a00  );
                    }
                  }
                }
              }
            } // end bounding box
          }  else if(as<std::string>(output)=="points"){
            // iterate over output points 
            for(int j=0; j<nG; j++) {
              ///////////////////////////////////////
              // calculate barycentric coordinates:
              double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[j] - yT[2])) /
                ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
              double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[j] - yT[2])) /
                ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
              double c = 1 - a - b;
              // check if inside triangle,
              // handle only yet untouched grid points
              // if(R_IsNA(z(j,0))){
              if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
                if(linear)
                  z(j,0)=a*zT[0]+b*zT[1]+c*zT[2];
                else{
                  // affine transformation
                  VectorXd xy(3), uv(3);
                  // use homogeneous coordinates:
                  xy << x[j], y[j], 1.0;
                  double uj,vj;
                    
                  uv = A*xy;
                  uj=uv[0];vj=uv[1];
                    
                  // no baryweight: only single term in summation
                  if(!baryweight)
                    a=1.0;
                    
                  if(i_barycycle==0){
                    z(j,0)=a*(  vj*(a01+uj*(uj*(a21+uj*(a41*uj+a31))+a11)+vj*(vj*(a03+uj*(a23*uj+a13)+vj*(a05*vj+a14*uj+a04))+uj*(a12+uj*(a32*uj+a22))+a02))+uj*(a10+uj*(uj*(a30+uj*(a50*uj+a40))+a20))+a00  );
                  } else {
                    z(j,0)+=a*(  vj*(a01+uj*(uj*(a21+uj*(a41*uj+a31))+a11)+vj*(vj*(a03+uj*(a23*uj+a13)+vj*(a05*vj+a14*uj+a04))+uj*(a12+uj*(a32*uj+a22))+a02))+uj*(a10+uj*(uj*(a30+uj*(a50*uj+a40))+a20))+a00  );
                  }
                }
              }
              ///////////////////////////////////////
            }
          } 
        } // i_barycycle loop
      } else { // !linear
      //################################### WORK 2 END
      iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
      xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
      yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];
      zT[0]=zD[tXY.i1[i]]; zT[1]=zD[tXY.i2[i]]; zT[2]=zD[tXY.i3[i]];


      if(as<std::string>(output)=="grid"){
	// iterate only over grid points (j,k) inside bounding box of triangle i
	for(int j=jTsw[i]; j<jTne[i]; j++) {
	  for(int k=kTsw[i]; k<kTne[i]; k++) {
	    // calculate barycentric coordinates:
	    double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[k] - yT[2])) /
	      ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	    double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[k] - yT[2])) /
	      ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	    double c = 1 - a - b;
	    // check if inside triangle, handle only yet untouched grid points
	    //if(R_IsNA(z(j,k))){
	    if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
	      z(j,k)=a*zT[0]+b*zT[1]+c*zT[2];
	    }
	  }
	}
      } else if(as<std::string>(output)=="points"){
	// iterate over output points
	for(int j=0; j<nG; j++) {
	  // calculate barycentric coordinates:
	  double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[j] - yT[2])) /
	    ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	  double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[j] - yT[2])) /
	    ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	  double c = 1 - a - b;
	  // check if inside triangle, handle only yet untouched grid points
	  //if(R_IsNA(z(j,k))){
	  if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
	    z(j,0)=a*zT[0]+b*zT[1]+c*zT[2];
	  }
	}
      }
            }// !linear

    } // triangle


    ret=List::create(_("x")=x, _("y")=y, _("z")=z);

        return ret;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached

}


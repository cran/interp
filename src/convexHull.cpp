// adopted from 
// https://algoteka.com/samples/35/graham-scan-convex-hull-algorithm-c-plus-plus-o%2528n-log-n%2529-readable-solution


/*
  LICENSE: https://algoteka.com/code-license

Copyright (c) 2022 Algoteka OÜ

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include<vector>
#include<algorithm>
#include<tuple>
#include<cmath>
#include <Rcpp.h>
using namespace Rcpp;

struct Point2D {
    double x;
    double y;
    
    Point2D (double cx, double cy){
      x=cx;
      y=cy;
    } 
    
    Point2D operator-(Point2D r) {
        return {x - r.x, y - r.y};
    }
    double operator*(Point2D r) {
        return x * r.x + y * r.y;
    }
    Point2D rotate90() {  // Rotate 90 degrees counter-clockwise
        return {-y, x};
    }
    double manhattan_length() {
        return fabs(x) + fabs(y);
    }
    bool operator==(Point2D r) {
        return x == r.x && y == r.y;
    }
    bool operator!=(Point2D r) {
        return x != r.x || y != r.y;
    }
};


std::vector<Point2D> graham_scan(std::vector<Point2D> points) {
    Point2D first_point = *std::min_element(points.begin(), points.end(), [](Point2D &left, Point2D &right) {
        return std::make_tuple(left.y, left.x) < std::make_tuple(right.y, right.x);
    });  // Find the lowest and leftmost point
    
    std::sort(points.begin(), points.end(), [&](Point2D &left, Point2D &right) {
        if(left == first_point) {
            return right != first_point;
        } else if (right == first_point) {
            return false;
        }
        double dir = (left-first_point).rotate90() * (right-first_point);
        if(dir == 0) {  // If the points are on a line with first point, sort by distance (manhattan is equivalent here)
            return (left-first_point).manhattan_length() < (right-first_point).manhattan_length();
        }
        return dir > 0;
        // Alternative approach, closer to common algorithm formulation but inferior:
        // return atan2(left.y - first_point.y, left.x - first_point.x) < atan2(right.y - first_point.y, right.x - first_point.x);
    });  // Sort the points by angle to the chosen first point
    
    std::vector<Point2D> result;
    for(auto pt : points) {
        // For as long as the last 3 points cause the hull to be non-convex, discard the middle one
        while (result.size() >= 2 &&
               (result[result.size()-1] - result[result.size()-2]).rotate90() * (pt - result[result.size()-1]) <= 0) {
            result.pop_back();
        }
        result.push_back(pt);
    }
    return result;
}

// [[Rcpp::export]]
List ConvexHull(NumericVector x, NumericVector y){
  
  int nx=x.size();
  int ny=y.size();
  List ret;
  std::vector<Point2D> pts;
  
    if(nx!=ny)
    ::Rf_error("ConvexHull: length of x and y dont match (%f!=%f)!",nx,ny);
    //Rcout << "prep" << std::endl;
    
    std::vector<double> vx=Rcpp::as<std::vector<double> >(x);
    std::vector<double> vy=Rcpp::as<std::vector<double> >(y);
    for(int i=0; i<nx; i++){
      pts.push_back(Point2D(vx[i],vy[i]));
    }
    //Rcout << "start" << std::endl;
    std::vector<Point2D> hull = graham_scan(pts);
    //Rcout << "extract" << std::endl;
    
    NumericVector hx(hull.size());
    NumericVector hy(hull.size());
    for(int i=0; i<hull.size();i++){
      hx[i]=hull[i].x;
      hy[i]=hull[i].y;
      //Rcout << i<< ": " << hx[i] << ", " <<hy[i]<<std::endl;
    }
    //Rcout << "ready:" << hx.size() << std::endl;
    
    return List::create(_("x")=hx, _("y")=hy);
}


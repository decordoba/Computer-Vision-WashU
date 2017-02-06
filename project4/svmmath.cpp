/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

#include "../Eigen/Core"
#include "MinEig.h"

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
  //Check if there are at least 2 lines
  if (lines.size() < 2)
  {
    fprintf(stderr, "Not enough lines to compute the best fit.");
    abort();
  }

  // To accumulate stuff
  typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

  // Transformation for numerical stability

  // Note: iterate through the lines list as follows:
  //		for (iter = lines.begin(); iter != lines.end(); iter++) {
  //			...iter is the pointer to the current line...
  //		}
  // Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
  //
  /******** BEGIN TODO ********/
  SVMPoint bestfit;
  Matrix3 M = Matrix3::Zero(3, 3);
  for (list<SVMLine>::const_iterator iter = lines.begin(); iter != lines.end(); iter++) {
    //Compute homogenous coordinate vector representing the line
    double a = iter->pnt1->v - iter->pnt2->v;
    double b = iter->pnt2->u - iter->pnt1->u;
    double c = iter->pnt1->u * iter->pnt2->v - iter->pnt1->v * iter->pnt2->u;
    //Form the 3x3 "second moment" matrix M
    M(0, 0) += a*a;
    M(0, 1) += a*b;
    M(0, 2) += a*c;
    M(1, 0) += a*b;
    M(1, 1) += b*b;
    M(1, 2) += b*c;
    M(2, 0) += a*c;
    M(2, 1) += b*c;
    M(2, 2) += c*c;
  }
  double eval, evec[3];
  //Perform an eigendecomposition of M
  MinEig(M, eval, evec);
  //Get the vanishing point
  bestfit.u = evec[0] / evec[2];
  bestfit.v = evec[1] / evec[2];
  bestfit.w = 1; //evec[2] / evec[2];
  /******** END TODO ********/
  return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		a plane coordinate of your choise.
//              See the pdf titled "Homography from Polygon in R^3 to Image Plane",
//              whose link can be found from the project page.
//
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
  int numPoints = points.size();
  if (numPoints < 3) {
    fprintf(stderr, "Not enough points to form a plane.");
    abort();
  }

  /******** BEGIN TODO ********/
//  for (int i = 0; i < numPoints; i++) {
//    printf("%d: (%f , %f, %f)  -  (%f , %f, %f, %f)\n", i, points.at(i).u, points.at(i).v, points.at(i).w, points.at(i).X, points.at(i).Y, points.at(i).Z, points.at(i).W);
//  }
//  printf("onverting to plane coordinate! zzzzzzzzzzzzzzzzzzzzzzzzzzz");

//  //Get first point (r)
//  Vec3d r = Vec3d(points.at(0).X, points.at(0).Y, points.at(0).Z);
//  //Get ex = (p-r)/nomr(p-r)
//  Vec3d ex = Vec3d(points.at(1).X, points.at(1).Y, points.at(1).Z) - r;
//  uScale = ex.length();
//  ex /= uScale;
//  //Get ey = (q-r)-(ex*dot_prod(q-r, ex))
//  Vec3d ey = Vec3d(points.at(2).X, points.at(2).Y, points.at(2).Z) - r;
//  double factor = ey[0] * ex[0] + ey[1] * ex[1] + ey[2] * ex[2];
//  ey -= (ex * factor);
//  vScale = ey.length();
//  ey /= vScale;

//  //Calculate every point's 2D coordinates from r, ex and ey
//  basisPts.clear();
//  for (int i = 0; i < numPoints; i++) {
//    Vec3d tmp = Vec3d(points.at(i).X, points.at(i).Y, points.at(i).Z) - r;
//    Vec3d pt = Vec3d(tmp[0] * ex[0] + tmp[1] * ex[1] + tmp[2] * ex[2], ey[0] * tmp[0] + ey[1] * tmp[1] + ey[2] * tmp[2], 1.0);
//    basisPts.push_back(pt);
//  }


  double best_prod = 0;
  Vec4d p = Vec4d(points[1].X, points[1].Y, points[1].Z, points[1].W);
  Vec4d q = Vec4d(points[2].X, points[2].Y, points[2].Z, points[2].W);
  Vec4d r = Vec4d(points[0].X, points[0].Y, points[0].Z, points[0].W);
  Vec4d pr = p - r;
  int pr_norm = sqrt(pr[0] * pr[0] + pr[1] * pr[1] + pr[2] * pr[2] + pr[3] * pr[3]);
  for(int i = 2; i < numPoints; i++){
    Vec4d qr = Vec4d(points[i].X, points[i].Y, points[i].Z, points[i].W) - r;
    int qr_norm = sqrt(qr[0]*qr[0] + qr[1]*qr[1] + qr[2]*qr[2] + qr[3]*qr[3]);
    double dot_prod = (pr[0]*qr[0] + pr[1]*qr[1] + pr[2]*qr[2] + pr[3]*qr[3]) / (qr_norm * pr_norm);
    if(dot_prod > 1){
      dot_prod = 1;
    }
    if(acos(dot_prod) > best_prod && acos(dot_prod) <= 3.1416/2){
      best_prod = acos(dot_prod);
      q = Vec4d(points[i].X, points[i].Y, points[i].Z, points[i].W);
    }
  }

  Vec4d qr = q - r;
  double ex_norm = pr.length();
  Vec4d ex = Vec4d(pr[0], pr[1], pr[2], pr[3]) / ex_norm;
  double dot_prod = ex[0] * qr[0] + ex[1] * qr[1] + ex[2] * qr[2] + ex[3] * qr[3];
  Vec4d s = Vec4d(ex[0],ex[1],ex[2],ex[3]) * dot_prod;
  Vec4d t = qr - s;
  double t_norm = t.length();
  Vec4d ey = Vec4d(t[0], t[1], t[2], t[3]) / t_norm;

  double min_u = std::numeric_limits<double>::max();
  double min_v = std::numeric_limits<double>::max();
  double max_u = std::numeric_limits<double>::min();
  double max_v = std::numeric_limits<double>::min();

  for (int i = 0; i < numPoints; i++) {
    Vec4d ar = Vec4d(points.at(i).X, points.at(i).Y, points.at(i).Z, points[i].W) - r;
    Vec3d pt = Vec3d(ar[0] * ex[0] + ar[1] * ex[1] + ar[2] * ex[2] + ar[3] * ex[3], ar[0] * ey[0] + ar[1] * ey[1] + ar[2] * ey[2] + ar[3] * ey[3], 1);
    basisPts.push_back(pt);
    if(min_u>basisPts[i][0]){
        min_u=basisPts[i][0];
    }
    if(max_u<basisPts[i][0]){
        max_u=basisPts[i][0];
    }
    if(min_v>basisPts[i][1]){
        min_v=basisPts[i][1];
    }
    if(max_u<basisPts[i][1]){
        max_v=basisPts[i][1];
    }
  }
  uScale=max_u-min_u;
  vScale=max_v-min_v;
  /******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see the pdf titled
//              "Homography from Polygon in R^3 to Image Plane", whose link can be found from
//              the project page.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
  int i;
  int numPoints = (int) points.size();

//  double uScaleD, vScaleD;
//  vector<Vec3d> basisPtsD;
//  vector<SVMPoint> pointsD;
//  SVMPoint p1;
//  p1.u = 200.5;
//  p1.v = 50;
//  p1.w = 1;
//  p1.X = 0;
//  p1.Y = 0;
//  p1.Z = 0;
//  p1.W = 1;
//  SVMPoint p2;
//  p2.u = 150.5;
//  p2.v = 75;
//  p2.w = 1;
//  p2.X = 0;
//  p2.Y = 1;
//  p2.Z = 0;
//  p2.W = 1;
//  SVMPoint p3;
//  p3.u = 250.5;
//  p3.v = 75;
//  p3.w = 1;
//  p3.X = 1;
//  p3.Y = 0;
//  p3.Z = 0;
//  p3.W = 1;
//  SVMPoint p4;
//  p4.u = 100.5;
//  p4.v = 100;
//  p4.w = 1;
//  p4.X = 0;
//  p4.Y = 1000;
//  p4.Z = 0;
//  p4.W = 1;
//  pointsD.push_back(p2);
//  pointsD.push_back(p1);
//  pointsD.push_back(p3);
//  pointsD.push_back(p4);

//  ConvertToPlaneCoordinate(pointsD, basisPtsD, uScaleD, vScaleD);
//  for (int i = 0; i < 4; i++) {
//    printf("%d --> (%f , %f, %f)\n", i, basisPtsD.at(i)[0], basisPtsD.at(i)[1], basisPtsD.at(i)[2]);
//  }


//  for (int i = 0; i < numPoints; i++) {
//    printf("%d: (%f , %f, %f)  -  (%f , %f, %f, %f)\n", i, points.at(i).u, points.at(i).v, points.at(i).w, points.at(i).X, points.at(i).Y, points.at(i).Z, points.at(i).W);
//  }
//  printf("ref plane: %d\n", isRefPlane);

  assert( numPoints >= 4 );

  basisPts.clear();
  if (isRefPlane) // reference plane
    {
      for (i=0; i < numPoints; i++)
        {
          Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
          basisPts.push_back(tmp);
          printf("%d: (%f , %f, %f)\n", i, tmp[0], tmp[1], tmp[2]);
        }
    }
  else // arbitrary polygon
    {
      double uScale, vScale; // unused in this function
      ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

  // A: 2n x 9 matrix where n is the number of points on the plane
  //    as discussed in lecture
  int numRows = 2 * numPoints;
  const int numCols = 9;

  typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
  MatrixType A = MatrixType::Zero(numRows, numCols);

  /******** BEGIN TODO ********/
  //Calculate umin, umax, vmin, vmax
//  double umax = std::numeric_limits<double>::min();
//  double umin = std::numeric_limits<double>::max();
//  double vmax = std::numeric_limits<double>::min();
//  double vmin = std::numeric_limits<double>::max();
//  for (i = 0; i < numPoints; i++) {
//    if (basisPts.at(i)[0] > umax) umax = basisPts.at(i)[0];
//    if (basisPts.at(i)[0] < umin) umin = basisPts.at(i)[0];
//    if (basisPts.at(i)[1] > vmax) vmax = basisPts.at(i)[1];
//    if (basisPts.at(i)[1] < vmin) vmin = basisPts.at(i)[1];
//  }
//  double udiff = umax-umin;
//  double vdiff = vmax-vmin;
//  printf("min max diff: %f, %f, %f, %f, %f, %f\n", umin, umax, vmin, vmax, udiff, vdiff);

//  for (i = 0; i < numPoints; i++) {
//    printf("%d: (%f, %f, %f)  -  (%f, %f, %f, %f)\n", i, points.at(i).u, points.at(i).v, points.at(i).w, points.at(i).X, points.at(i).Y, points.at(i).Z, points.at(i).W);
//  }
//  for (i = 0; i < numPoints; i++) {
//    printf("%d: (%f, %f, %f)\n", i, basisPts.at(i)[0], basisPts.at(i)[1], basisPts.at(i)[2]);
//  }


  //Fill A
  for (i = 0; i < numPoints; i++) {
//    double x = - (basisPts.at(i)[0] - umin) / udiff;
//    double y = - (basisPts.at(i)[1] - vmin) / vdiff;
//    double x = - basisPts.at(i)[0];
//    double y = - basisPts.at(i)[1];
//    printf("%d: %f, %f\n", i, x, y);
//    A(i*2, 0) = points.at(i).v;
//    A(i*2, 1) = points.at(i).u;
//    A(i*2, 2) = points.at(i).w;
//    A(i*2, 6) = points.at(i).v * x;
//    A(i*2, 7) = points.at(i).u * x;
//    A(i*2, 8) = points.at(i).w * x;
//    A(i*2+1, 3) = points.at(i).v;
//    A(i*2+1, 4) = points.at(i).u;
//    A(i*2+1, 5) = points.at(i).w;
//    A(i*2+1, 6) = points.at(i).v * y;
//    A(i*2+1, 7) = points.at(i).u * y;
//    A(i*2+1, 8) = points.at(i).w * y;

    A(i*2, 0) = basisPts.at(i)[0];
    A(i*2, 1) = basisPts.at(i)[1];
    A(i*2, 2) = 1;
    A(i*2, 6) = - points.at(i).u * basisPts.at(i)[0];
    A(i*2, 7) = - points.at(i).u * basisPts.at(i)[1];
    A(i*2, 8) = - points.at(i).u;
    A(i*2+1, 3) = basisPts.at(i)[0];
    A(i*2+1, 4) = basisPts.at(i)[1];
    A(i*2+1, 5) = 1;
    A(i*2+1, 6) = - points.at(i).v * basisPts.at(i)[0];
    A(i*2+1, 7) = - points.at(i).v * basisPts.at(i)[1];
    A(i*2+1, 8) = - points.at(i).v;
  }

  for (i = 0; i < numRows; i++) {
      printf("%d: (%f, %f, %f, %f, %f, %f, %f, %f, %f)\n", i, A(i,0), A(i,1), A(i,2), A(i,3), A(i,4), A(i,5), A(i,6), A(i,7), A(i,8));
  }

  double eval, h[9];
  MinEig(A, eval, h);

  H[0][0] = h[0];
  H[0][1] = h[1];
  H[0][2] = h[2];

  H[1][0] = h[3];
  H[1][1] = h[4];
  H[1][2] = h[5];

  H[2][0] = h[6];
  H[2][1] = h[7];
  H[2][2] = h[8];
  /******** END TODO ********/

  // compute inverse of H
  if (H.Determinant() == 0)
    fl_alert("Computed homography matrix is uninvertible \n");
  else
    Hinv = H.Inverse();

  int ii;
  printf("\nH=[\n");
  for (ii=0; ii<3; ii++)
    printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
  printf("]\nHinv=[\n");

  for (ii=0; ii<3; ii++)
    printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);



  printf("]\n\n");
}

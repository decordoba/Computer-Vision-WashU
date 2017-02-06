/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
    if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  if( refPointOffPlane == NULL )
    {
      fl_alert("Need to specify the reference height first.");
      return;
    }

  printf("Ref (%f , %f, %f)  -  (%f , %f, %f, %f)\n", refPointOffPlane->u, refPointOffPlane->v, refPointOffPlane->w, refPointOffPlane->X, refPointOffPlane->Y, refPointOffPlane->Z, refPointOffPlane->W);
  
  /******** BEGIN TODO ********/
  
  // See the lecture note on measuring heights
  // using a known point directly below the new point.

  //Normalize everything
  Vec3d refP = Vec3d(refPointOffPlane->u/refPointOffPlane->w,refPointOffPlane->v/refPointOffPlane->w, 1);
  Vec3d vz = Vec3d(zVanish.u/zVanish.w, zVanish.v/zVanish.w, 1);
  Vec3d newP = Vec3d(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);
  Vec3d knownP = Vec3d(knownPoint.u/knownPoint.w, knownPoint.v/knownPoint.w, 1);
  //Calculate part equation
  double tb = (knownP - newP) * (knownP - newP);
  double rb = (knownP - refP) * (knownP - refP);
  double vzr = (vz - refP) * (vz - refP);
  double vzt = (vz - newP) * (vz - newP);
  //Save new point
  newPoint.Z = referenceHeight * sqrt((tb * vzr) / (rb * vzt));
  newPoint.X = knownPoint.X;
  newPoint.Y = knownPoint.Y;
  newPoint.W = 1;
  /******** END TODO ********/
 
  newPoint.known(true);
 
  printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
  redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  /******** BEGIN TODO ********/
  Vec3d b1 = Vec3d(newPoint.u / newPoint.w, newPoint.v / newPoint.w, 1);
  if (knownPoint.Z != 0){
    //New and known points
    Vec3d newP = Vec3d(newPoint.u, newPoint.v, newPoint.w);
    Vec3d knownP = Vec3d(knownPoint.u, knownPoint.v, knownPoint.w);
    //Vanishing points
    Vec3d vx = Vec3d(xVanish.u, xVanish.v, xVanish.w);
    Vec3d vy = Vec3d(yVanish.u, yVanish.v, yVanish.w);
    Vec3d vz = Vec3d(zVanish.u, zVanish.v, zVanish.w);
    //Horizon
    Vec3d horizon = cross(vx, vy);
    //Vertical line that goes through new point
    Vec3d vertical = cross(newP, vz);
    //Vanishing point for known and new points
    Vec3d v = cross(cross(newP, knownP), horizon);
    if (v * v == 0) {
      //Point of intersection horizon - vertical
      b1 = cross(horizon, vertical);
    } else {
      Mat3d H_mat = Mat3d(H[0][0],H[0][1],H[0][2],
                          H[1][0],H[1][1],H[1][2],
                          H[2][0],H[2][1],H[2][2]);
      //Image coordinates of the known point from homography
      Vec3d b0 = H_mat * Vec3d(knownPoint.X, knownPoint.Y, 1);
      b1 = cross(cross(b0, v), vertical);
    }
    //Normalize
    b1[0] = b1[0] / b1[2];
    b1[1] = b1[1] / b1[2];
    b1[2] = 1;
  }
  Mat3d Hinv_mat = Mat3d(Hinv[0][0], Hinv[0][1], Hinv[0][2],
                         Hinv[1][0], Hinv[1][1], Hinv[1][2],
                         Hinv[2][0], Hinv[2][1], Hinv[2][2]);
  Vec3d updatedP = Hinv_mat * b1;
  //Save point
  newPoint.X = updatedP[0] / updatedP[2];
  newPoint.Y = updatedP[1] / updatedP[2];
  newPoint.Z = knownPoint.Z;
  newPoint.W = 1;
  /******** END TODO ********/
 
  newPoint.known(true);

  printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
  redraw();
}


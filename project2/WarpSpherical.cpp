///////////////////////////////////////////////////////////////////////////
//
// NAME
//  WarpSpherical.h -- warp a flat (perspective) image into spherical
//      coordinates and/or undo radial lens distortion
//
// SEE ALSO
//  WarpSpherical.h   longer description
//
//  R. Szeliski and H.-Y. Shum.
//  Creating full view panoramic image mosaics and texture-mapped models.
//  Computer Graphics (SIGGRAPH'97), pages 251-258, August 1997.
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE455 Spring 2008)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "WarpSpherical.h"
#include <math.h>

/******************* TO DO *********************
 * warpSphericalField:
 *	INPUT:
 *		srcSh: shape (width, height, number of channels) of source image
 *		dstSh: shape of destination image
 *		f: focal length in pixels (provided on the project web page, or measured by yourself)
 *		k1, k2: radial distortion parameters (ditto)
 *		r: rotation matrix
 *	OUTPUT:
 *		Return an image containing (u,v) coordinates for mapping pixels from
 *		spherical coordinates to planar image coordinates and applying
 *		radial distortion.
 *			Note that this is inverse warping, i.e., this routine will be
 *		actually used to map from planar coordinates with radial distortion
 *		to spherical coordinates without radial distortion.
 *
 */
CFloatImage WarpSphericalField(CShape srcSh, CShape dstSh, float f,
                                 float k1, float k2, const CTransform3x3 &r)
{
    // Set up the pixel coordinate image
    dstSh.nBands = 2;
    CFloatImage uvImg(dstSh);   // (u,v) coordinates

    // Fill in the values
    for (int y = 0; y < dstSh.height; y++)
    {
        float *uv = &uvImg.Pixel(0, y, 0);
        for (int x = 0; x < dstSh.width; x++, uv += 2)
        {
			// (x,y) is the spherical image coordinates. 
            // (xf,yf) is the spherical coordinates, e.g., xf is the angle theta
			// and yf is the angle phi

            float xf = (x - 0.5f*dstSh.width ) / f;
            float yf = (y - 0.5f*dstSh.height) / f;

			// (xt,yt,zt) are intermediate coordinates to which you can
			// apply the spherical correction and radial distortion
            float xt, yt, zt;
			CVector3 p;
 
			// *** BEGIN TODO ***
			// add code to apply the spherical correction, i.e.,
			// compute the Euclidean coordinates, rotate according to
			// r, and project the point to the z=1 plane at
			// (xt/zt,yt/zt,1), then distort with radial distortion
            // coefficients k1 and k2

            // Compute Euclidean coordinates
            p[0] = sin(xf) * cos(yf);
            p[1] = sin(yf);
            p[2] = cos(xf) * cos(yf);

            // Apply rotation
            // xt = r[0][0]*p[0] + r[0][1]*p[1] + r[0][2]*p[2];
            // yt = r[1][0]*p[0] + r[1][1]*p[1] + r[1][2]*p[2];
            // zt = r[2][0]*p[0] + r[2][1]*p[1] + r[2][2]*p[2];
            p = r*p;

            // Project to z=1 plane
            xt = p[0]/p[2];
            yt = p[1]/p[2];
            zt = p[2]/p[2];

            // Distort with radial distortion
            float r2 = pow(xt, 2) + pow(yt, 2);
            float r4 = pow(r2, 2);
            xt = xt * (1 + k1*r2 + k2*r4);
            yt = yt * (1 + k1*r2 + k2*r4);
			// *** END TODO ***

            // Convert back to regular pixel coordinates and store
            float xn = 0.5f*srcSh.width  + xt*f;
            float yn = 0.5f*srcSh.height + yt*f;
            uv[0] = xn;
            uv[1] = yn;
        }
    }
    return uvImg;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BlendImages.cpp -- blend together a set of overlapping images
//
// DESCRIPTION
//  This routine takes a collection of images aligned more or less horizontally
//  and stitches together a mosaic.
//
//  The images can be blended together any way you like, but I would recommend
//  using a soft halfway blend of the kind Steve presented in the first lecture.
//
//  Once you have blended the images together, you should crop the resulting
//  mosaic at the halfway points of the first and last image.  You should also
//  take out any accumulated vertical drift using an affine warp.
//  Lucas-Kanade Taylor series expansion of the registration error.
//
// SEE ALSO
//  BlendImages.h       longer description of parameters
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE455 Winter 2003)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "BlendImages.h"
#include <float.h>
#include <math.h>

#include <stdio.h>

#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

/* Return the closest integer to x, rounding up */
static int iround(double x) {
    if (x < 0.0) {
	return (int) (x - 0.5);
    } else {
	return (int) (x + 0.5);
    }
}

void ImageBoundingBox(CImage &image, CTransform3x3 &M, 
		      int &min_x, int &min_y, int &max_x, int &max_y)
{
    float min_xf = FLT_MAX, min_yf = FLT_MAX;
    float max_xf = 0.0, max_yf = 0.0;

    CVector3 corners[4];

    int width = image.Shape().width;
    int height = image.Shape().height;
	
    corners[0][0] = 0.0;
    corners[0][1] = 0.0;
    corners[0][2] = 1.0;
	
    corners[1][0] = width - 1;
    corners[1][1] = 0.0;
    corners[1][2] = 1.0;

    corners[2][0] = 0.0;
    corners[2][1] = height - 1;
    corners[2][2] = 1.0;

    corners[3][0] = width - 1;
    corners[3][1] = height - 1;
    corners[3][2] = 1.0;

    corners[0] = M * corners[0];
    corners[1] = M * corners[1];
    corners[2] = M * corners[2];
    corners[3] = M * corners[3];

    corners[0][0] /= corners[0][2];
    corners[0][1] /= corners[0][2];

    corners[1][0] /= corners[0][2];
    corners[1][1] /= corners[0][2];

    corners[2][0] /= corners[0][2];
    corners[2][1] /= corners[0][2];

    corners[3][0] /= corners[0][2];
    corners[3][1] /= corners[0][2];
	

    min_xf = (float) MIN(min_xf, corners[0][0]);
    min_xf = (float) MIN(min_xf, corners[1][0]);
    min_xf = (float) MIN(min_xf, corners[2][0]);
    min_xf = (float) MIN(min_xf, corners[3][0]);

    min_yf = (float) MIN(min_yf, corners[0][1]);
    min_yf = (float) MIN(min_yf, corners[1][1]);
    min_yf = (float) MIN(min_yf, corners[2][1]);
    min_yf = (float) MIN(min_yf, corners[3][1]);
	
    max_xf = (float) MAX(max_xf, corners[0][0]);
    max_xf = (float) MAX(max_xf, corners[1][0]);
    max_xf = (float) MAX(max_xf, corners[2][0]);
    max_xf = (float) MAX(max_xf, corners[3][0]);
	
    max_yf = (float) MAX(max_yf, corners[0][1]);
    max_yf = (float) MAX(max_yf, corners[1][1]);
    max_yf = (float) MAX(max_yf, corners[2][1]);
    max_yf = (float) MAX(max_yf, corners[3][1]);    

    min_x = (int) floor(min_xf);
    min_y = (int) floor(min_yf);
    max_x = (int) ceil(max_xf);
    max_y = (int) ceil(max_yf);
}


/******************* TO DO *********************
 * AccumulateBlend:
 *	INPUT:
 *		img: a new image to be added to acc
 *		acc: portion of the accumulated image where img is to be added
 *		M: translation matrix for calculating a bounding box
 *		blendWidth: width of the blending function (horizontal hat function;
 *	    try other blending functions for extra credit)
 *	OUTPUT:
 *		add a weighted copy of img to the subimage specified in acc
 *		the first 3 band of acc records the weighted sum of pixel colors
 *		the fourth band of acc records the sum of weight
 */
static void AccumulateBlend(CByteImage& img, CFloatImage& acc, CTransform3x3 M, float blendWidth)
{
    /* Compute the bounding box of the image of the image */
    int bb_min_x, bb_min_y, bb_max_x, bb_max_y;
    ImageBoundingBox(img, M, bb_min_x, bb_min_y, bb_max_x, bb_max_y);

    // Used to show blending sections
    bool see_blend = false;
    // Used to activate optimal blending
    int optimal_blend = 1; // 0:simple, 1:optimal, 2:average
    // Used to see images for debugging and pause the program every iteration
    bool debug_mode = false;
    // Correct exposure differences
    bool color_correct = false;

    CTransform3x3 Minv = M.Inverse();

    // Find start and end between intersection with previous pictures
    // Also correct colors
    int intersection_start = acc.Shape().width;
    int intersection_end = 0;
    int R1avg = 0, R2avg = 0, G1avg = 0, G2avg = 0, B1avg = 0, B2avg = 0;
    int numRGB = 0;
    for (int x = bb_min_x; x < bb_max_x; x++) {
        for (int y = bb_min_y; y <= bb_max_y; y++) {
            if (acc.Pixel(x, y, 3) > 0) {
                if (intersection_end < x)
                    intersection_end = x;
                if (intersection_start > x)
                    intersection_start = x;
            }

            CVector3 p_dest, p_src;
            p_dest[0] = x;
            p_dest[1] = y;
            p_dest[2] = 1.0;
            p_src = Minv * p_dest;

            float x_src = (float) (p_src[0] / p_src[2]);
            float y_src = (float) (p_src[1] / p_src[2]);
            if (x_src < 0.0 || x_src >= img.Shape().width - 1 || y_src < 0.0 || y_src >= img.Shape().height - 1)
                continue;
            int xf = (int) floor(x_src);
            int yf = (int) floor(y_src);
            int xc = xf + 1;
            int yc = yf + 1;
            /* Skip black pixels */
            if (img.Pixel(xf, yf, 0) == 0x0 && img.Pixel(xf, yf, 1) == 0x0 && img.Pixel(xf, yf, 2) == 0x0)
                continue;
            if (img.Pixel(xc, yf, 0) == 0x0 && img.Pixel(xc, yf, 1) == 0x0 && img.Pixel(xc, yf, 2) == 0x0)
                continue;
            if (img.Pixel(xf, yc, 0) == 0x0 && img.Pixel(xf, yc, 1) == 0x0 && img.Pixel(xf, yc, 2) == 0x0)
                continue;
            if (img.Pixel(xc, yc, 0) == 0x0 && img.Pixel(xc, yc, 1) == 0x0 && img.Pixel(xc, yc, 2) == 0x0)
                continue;
            if (acc.Pixel(x, y, 3) > 0) {
                // If pixel exists in acc and img, store RGB
                R1avg += acc.Pixel(x, y, 0);
                G1avg += acc.Pixel(x, y, 1);
                B1avg += acc.Pixel(x, y, 2);
                numRGB++;
                R2avg += img.PixelLerp(x_src, y_src, 0);
                G2avg += img.PixelLerp(x_src, y_src, 1);
                B2avg += img.PixelLerp(x_src, y_src, 2);
            }
        }
    }
    // Compute RGB difference between two images
    float Roffset = 0;
    float Goffset = 0;
    float Boffset = 0;
    if (numRGB > 0 && color_correct) {
        Roffset = (R2avg/((float) numRGB) - R1avg/((float) numRGB));
        Goffset = (G2avg/((float) numRGB) - G1avg/((float) numRGB));
        Boffset = (B2avg/((float) numRGB) - B1avg/((float) numRGB));
    }

    // Middle of both intersections
    int intersection_middle = (intersection_start + intersection_end) / 2;
    // Bool that says if images are pasted to the right or left
    bool paste_right = (intersection_start - bb_min_x) < (bb_max_x - intersection_end);

    for (int y = bb_min_y; y <= bb_max_y; y++) {
        for (int x = bb_min_x; x < bb_max_x; x++) {
            /* Check bounds in destination */
            if (x < 0 || x >= acc.Shape().width || y < 0 || y >= acc.Shape().height) {
                continue;
            }

            /* Compute source pixel and check bounds in source */
            CVector3 p_dest, p_src;
            p_dest[0] = x;
            p_dest[1] = y;
            p_dest[2] = 1.0;

            p_src = Minv * p_dest;

            float x_src = (float) (p_src[0] / p_src[2]);
            float y_src = (float) (p_src[1] / p_src[2]);

            if (x_src < 0.0 || x_src >= img.Shape().width - 1 ||
                y_src < 0.0 || y_src >= img.Shape().height - 1) {
                continue;
            }

            int xf = (int) floor(x_src);
            int yf = (int) floor(y_src);
            int xc = xf + 1;
            int yc = yf + 1;

            /* Skip black pixels */
            if (img.Pixel(xf, yf, 0) == 0x0 &&
                img.Pixel(xf, yf, 1) == 0x0 &&
                img.Pixel(xf, yf, 2) == 0x0) {
                continue;
            }

            if (img.Pixel(xc, yf, 0) == 0x0 &&
                img.Pixel(xc, yf, 1) == 0x0 &&
                img.Pixel(xc, yf, 2) == 0x0) {
                continue;
            }

            if (img.Pixel(xf, yc, 0) == 0x0 &&
                img.Pixel(xf, yc, 1) == 0x0 &&
                img.Pixel(xf, yc, 2) == 0x0) {
                continue;
            }

            if (img.Pixel(xc, yc, 0) == 0x0 &&
                img.Pixel(xc, yc, 1) == 0x0 &&
                img.Pixel(xc, yc, 2) == 0x0) {
                continue;
            }

            // Compute weight to perform blending
            double weight = 1.0;

            // *** BEGIN TODO ***
            // set weight properly

            // If we are looking at a non transparent pixel in acc (we are intersecting)
            if (acc.Pixel(x, y, 3) > 0)
            {
                if (paste_right) {
                    if (optimal_blend == 2) {
                        // Regular feathering
                        if (x - bb_min_x < blendWidth) {
                            weight = (x - bb_min_x) / blendWidth;
                        }
                    }
                    else {
                        // Optimal blend
                        if (x < intersection_middle - blendWidth/2) {
                            // From left of intersection_middle - blendWidth/2 we only show acc
                            weight = 0.0;
                        }
                        else if (x > intersection_middle + blendWidth/2) {
                            // From right of intersection_middle + blendWidth/2 we only show current img
                            weight = 1.0;
                        }
                        else {
                            // Else we blend with weights depending on distance from intersection_middle
                            weight = (x - intersection_middle + blendWidth/2) / blendWidth;
                        }
                    }
                }
                else {
                    if (optimal_blend == 2) {
                        // Regular feathering
                        if (bb_max_x - x < blendWidth) {
                            weight = (bb_max_x - x) / blendWidth;
                        }
                    }
                    else {
                        // Optimal blend
                        if (x < intersection_middle - blendWidth/2) {
                            // From left of intersection_middle - blendWidth/2 we only show current img
                            weight = 1.0;
                        }
                        else if (x > intersection_middle + blendWidth/2) {
                            // From right of intersection_middle + blendWidth/2 we only show acc
                            weight = 0.0;
                        }
                        else {
                            // Else we blend with weights depending on distance from intersection_middle
                            weight = -(x - intersection_middle - blendWidth/2) / blendWidth;
                        }
                    }
                }
            }

            if (optimal_blend > 0) {
                acc.Pixel(x, y, 0) = (float) ((1-weight) * acc.Pixel(x, y, 0)) + (float) (weight * MAX(0, MIN(255.0, img.PixelLerp(x_src, y_src, 0) - Roffset)));
                acc.Pixel(x, y, 1) = (float) ((1-weight) * acc.Pixel(x, y, 1)) + (float) (weight * MAX(0, MIN(255.0, img.PixelLerp(x_src, y_src, 1) - Goffset)));
                acc.Pixel(x, y, 2) = (float) ((1-weight) * acc.Pixel(x, y, 2)) + (float) (weight * MAX(0, MIN(255.0, img.PixelLerp(x_src, y_src, 2) - Boffset)));
                acc.Pixel(x, y, 3) = (float) 1;
            }
            else {
                // This is the method for simple blend
                weight = 1.0;
                acc.Pixel(x, y, 0) += (float) (weight * img.PixelLerp(x_src, y_src, 0));
                acc.Pixel(x, y, 1) += (float) (weight * img.PixelLerp(x_src, y_src, 1));
                acc.Pixel(x, y, 2) += (float) (weight * img.PixelLerp(x_src, y_src, 2));
                acc.Pixel(x, y, 3) += weight;
            }

            // Used to see the blend region between images and visually determine the optimal blendWidth
            if (see_blend) {
                if (x < intersection_middle+blendWidth/2 && x > intersection_middle-blendWidth/2) {
                    acc.Pixel(x, y, 0) = 255.0;
                    acc.Pixel(x, y, 1) = 0.0;
                    acc.Pixel(x, y, 2) = 0.0;
                    acc.Pixel(x, y, 1) = 1.0;
                }
                if (x >= intersection_middle-2 && x <= intersection_middle+2) {
                    acc.Pixel(x, y, 0) = 255.0;
                    acc.Pixel(x, y, 1) = 255.0;
                    acc.Pixel(x, y, 2) = 255.0;
                    acc.Pixel(x, y, 1) = 1.0;
                }
            }
            // *** END TODO ***
        }
    }

    // For debugging purposes
    if (debug_mode) {
        CByteImage imgtmp(acc.Shape());

        uchar min = 10000, max = 0;
        uchar min2 = 10000, max2 = 0;
        for (int y = 0; y < acc.Shape().height; y++) {
            for (int x = 0; x < acc.Shape().width; x++) {
                imgtmp.Pixel(x, y, 0) = (uchar) floor(acc.Pixel(x, y, 0) / acc.Pixel(x, y, 3));
                imgtmp.Pixel(x, y, 1) = (uchar) floor(acc.Pixel(x, y, 1) / acc.Pixel(x, y, 3));
                imgtmp.Pixel(x, y, 2) = (uchar) floor(acc.Pixel(x, y, 2) / acc.Pixel(x, y, 3));
                imgtmp.Pixel(x, y, 3) = (uchar) 255.0;
                if (min > imgtmp.Pixel(x, y, 0)) min = imgtmp.Pixel(x, y, 0);
                if (max < imgtmp.Pixel(x, y, 0)) max = imgtmp.Pixel(x, y, 0);
                if (min2 > imgtmp.Pixel(x, y, 3)) min2 = imgtmp.Pixel(x, y, 3);
                if (max2 < imgtmp.Pixel(x, y, 3)) max2 = imgtmp.Pixel(x, y, 3);
            }
        }

        for (int y = 0; y < img.Shape().height; y++) {
            for (int x = 0; x < img.Shape().width; x++) {
                if (min2 > img.Pixel(x, y, 3)) min2 = img.Pixel(x, y, 3);
                if (max2 < img.Pixel(x, y, 3)) max2 = img.Pixel(x, y, 3);
            }
        }

        WriteFile(img, "my_img.tga");
        printf("my_img.tga saved\n");

        WriteFile(imgtmp, "my_new_img.tga");
        printf("my_new_img.tga saved\n");

        printf("Press Enter to continue\n");
        getchar();
    }
}


/******************* TO DO *********************
 * NormalizeBlend:
 *	INPUT:
 *		acc: input image whose alpha channel (4th channel) contains
 *		     normalizing weight values
 *		img: where output image will be stored
 *	OUTPUT:
 *		normalize r,g,b values (first 3 channels) of acc and store it into img
 */
static void NormalizeBlend(CFloatImage& acc, CByteImage& img)
{
	// *** BEGIN TODO ***
    // fill in this routine..
    for (int y = 0; y < acc.Shape().height; y++) {
        for (int x = 0; x < acc.Shape().width; x++) {
            if (acc.Pixel(x, y, 3) < 1) {
                img.Pixel(x, y, 3) = (uchar) 0;
            } else {
                img.Pixel(x, y, 0) = (uchar) (acc.Pixel(x, y, 0) / acc.Pixel(x, y, 3));
                img.Pixel(x, y, 1) = (uchar) (acc.Pixel(x, y, 1) / acc.Pixel(x, y, 3));
                img.Pixel(x, y, 2) = (uchar) (acc.Pixel(x, y, 2) / acc.Pixel(x, y, 3));
                img.Pixel(x, y, 3) = (uchar) 255;
            }
        }
    }
    // *** END TODO ***
}


/******************* TO DO *********************
 * BlendImages:
 *	INPUT:
 *		ipv: list of input images and their relative positions in the mosaic
 *		blendWidth: width of the blending function
 *	OUTPUT:
 *		create & return final mosaic by blending all images
 *		and correcting for any vertical drift
 */
CByteImage BlendImages(CImagePositionV& ipv, float blendWidth)
{
    // Assume all the images are of the same shape (for now)
    CByteImage& img0 = ipv[0].img;
    CShape sh        = img0.Shape();
    int width        = sh.width;
    int height       = sh.height;
    int nBands       = sh.nBands;
    int dim[2]       = {width, height};

    // Compute the bounding box for the mosaic
    int n = ipv.size();
    float min_x = 0, min_y = 0;
    float max_x = 0, max_y = 0;

    CVector3 lastCorner[4];
    CVector3 firstCorner[4];

    for (int i = 0; i < n; i++)
    {
		
		CTransform3x3 &pos = ipv[i].position;

        CVector3 corners[4];

        corners[0][0] = 0.0;
        corners[0][1] = 0.0;
        corners[0][2] = 1.0;

        corners[1][0] = width - 1;
        corners[1][1] = 0.0;
        corners[1][2] = 1.0;

        corners[2][0] = 0.0;
        corners[2][1] = height - 1;
        corners[2][2] = 1.0;

        corners[3][0] = width - 1;
        corners[3][1] = height - 1;
        corners[3][2] = 1.0;

        corners[0] = pos * corners[0];
        corners[1] = pos * corners[1];
        corners[2] = pos * corners[2];
        corners[3] = pos * corners[3];

        corners[0][0] /= corners[0][2];
        corners[0][1] /= corners[0][2];

        corners[1][0] /= corners[0][2];
        corners[1][1] /= corners[0][2];

        corners[2][0] /= corners[0][2];
        corners[2][1] /= corners[0][2];

        corners[3][0] /= corners[0][2];
        corners[3][1] /= corners[0][2];
        
        // *** BEGIN TODO #1 ***
        // add some code here to update min_x, ..., max_y

        if (i == 0) {
            firstCorner[0][0] = corners[0][0];
            firstCorner[0][1] = corners[0][1];
            firstCorner[0][2] = corners[0][2];
            firstCorner[1][0] = corners[1][0];
            firstCorner[1][1] = corners[1][1];
            firstCorner[1][2] = corners[1][2];
            firstCorner[2][0] = corners[2][0];
            firstCorner[2][1] = corners[2][1];
            firstCorner[2][2] = corners[2][2];
            firstCorner[3][0] = corners[3][0];
            firstCorner[3][1] = corners[3][1];
            firstCorner[3][2] = corners[3][2];
        }
        if (i == n-1) {
            lastCorner[0][0] = corners[1][0];
            lastCorner[0][1] = corners[1][1];
            lastCorner[0][2] = corners[1][2];
            lastCorner[1][0] = corners[0][0];
            lastCorner[1][1] = corners[0][1];
            lastCorner[1][2] = corners[0][2];
            lastCorner[2][0] = corners[3][0];
            lastCorner[2][1] = corners[3][1];
            lastCorner[2][2] = corners[3][2];
            lastCorner[3][0] = corners[2][0];
            lastCorner[3][1] = corners[2][1];
            lastCorner[3][2] = corners[2][2];
        }

        min_x = (float) MIN(min_x, corners[0][0]);
        min_x = (float) MIN(min_x, corners[1][0]);
        min_x = (float) MIN(min_x, corners[2][0]);
        min_x = (float) MIN(min_x, corners[3][0]);

        min_y = (float) MIN(min_y, corners[0][1]);
        min_y = (float) MIN(min_y, corners[1][1]);
        min_y = (float) MIN(min_y, corners[2][1]);
        min_y = (float) MIN(min_y, corners[3][1]);

        max_x = (float) MAX(max_x, corners[0][0]);
        max_x = (float) MAX(max_x, corners[1][0]);
        max_x = (float) MAX(max_x, corners[2][0]);
        max_x = (float) MAX(max_x, corners[3][0]);

        max_y = (float) MAX(max_y, corners[0][1]);
        max_y = (float) MAX(max_y, corners[1][1]);
        max_y = (float) MAX(max_y, corners[2][1]);
        max_y = (float) MAX(max_y, corners[3][1]);
		// *** END TODO #1 ***
    }

    // Create a floating point accumulation image
    CShape mShape((int)(ceil(max_x) - floor(min_x)),
                  (int)(ceil(max_y) - floor(min_y)), nBands);
    CFloatImage accumulator(mShape);
    accumulator.ClearPixels();

    printf("Shape accumulator: %d x %d\n", mShape.width, mShape.height);

	double x_init, x_final;
    double y_init, y_final;

	// Add in all of the images
    for (int i = 0; i < n; i++) {
        
        CTransform3x3 &M = ipv[i].position;

        CTransform3x3 M_t = CTransform3x3::Translation(-min_x, -min_y) * M;

        CByteImage& img = ipv[i].img;

        // Perform the accumulation
        AccumulateBlend(img, accumulator, M_t, blendWidth);

        if (i == 0) {
            CVector3 p;
            p[0] = 0.5 * width;
            p[1] = 0.0;
            p[2] = 1.0;

            p = M_t * p;
            x_init = p[0];
            y_init = p[1];
        } else if (i == n - 1) {
            CVector3 p;
            p[0] = 0.5 * width;
            p[1] = 0.0;
            p[2] = 1.0;

            p = M_t * p;
            x_final = p[0];
            y_final = p[1];
        }
    }

    // Normalize the results
    CByteImage compImage(mShape);
    NormalizeBlend(accumulator, compImage);
    bool debug_comp = true;
    if (debug_comp) {
        printf("Debug image tmp_comp.tga created\n");
        WriteFile(compImage, "tmp_comp.tga");
    }

    // Allocate the final image shape

//    // Yosemite (blendWidth=150)
//    float margin_x = 0.02, margin_y = 0.11; // Margin to remove
//    float offset_x = 0.00, offset_y = -0.03; // Offsets for x and y

//    // Campus (blendWidth=200)
//    float margin_x = 0.005, margin_y = 0.14; // Margin to remove
//    float offset_x = 0.00, offset_y = 0.05; // Offsets for x and y

//    // Room (blendWidth=100)
//    float margin_x = 0.015, margin_y = 0.22; // Margin to remove
//    float offset_x = 0.00, offset_y = -0.065; // Offsets for x and y

    // No cropping
    float margin_x = 0.00, margin_y = 0.00; // Margin to remove
    float offset_x = 0.00, offset_y = 0.00; // Offsets for x and y

    CShape cShape(floor(mShape.width*(1.0-2.0*margin_x)), floor(mShape.height*(1.0-2.0*margin_y)), nBands);
    printf("Final size: %dx%d\n", cShape.width, cShape.height);
    CByteImage croppedImage(cShape);

    // Compute the affine deformation
    CTransform3x3 A;
    
	// *** BEGIN TODO #2 ***
    // fill in the right entries in A to trim the left edge and
    // to take out the vertical drift

    A[0][2] = (margin_x + offset_x) * mShape.width ;
    A[1][2] = (margin_y + offset_y) * mShape.height;

    // Different distance if panorama was made turning left or right
    float distCornersX = lastCorner[0][0] - firstCorner[0][0];
    if (abs(distCornersX) < abs(lastCorner[1][0] - firstCorner[1][0])) {
        distCornersX = lastCorner[1][0] - firstCorner[1][0];
    }
    A[1][0] = (lastCorner[0][1]-firstCorner[0][1]) / (distCornersX);

	// *** END TODO #2 ***

    // Warp and crop the composite
    WarpGlobal(compImage, croppedImage, A, eWarpInterpLinear);

    return croppedImage;
}

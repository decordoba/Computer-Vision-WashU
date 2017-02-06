#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.
	switch (featureType) {
	case 1:
		dummyComputeFeatures(image, features);
		break;
	case 2:
		ComputeHarrisFeatures(image, features);
		break;
	default:
		return false;
	}

	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i=0; i<features.size(); i++) {
		features[i].id = i+1;
	}

	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i=0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.
	
	printf("\nMatching features.......\n");

	switch (matchType) {
    case 1:
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
        ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
        return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
		d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		n++;
	}	

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9], vector<bool> &isMatch, double threshold, double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
		applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		if (d<=threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD=matches[i].score;
	}	
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i=0; i < (int)thresholds.size();i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp=0;
		int actualCorrect=0;
		int fp=0;
		int actualError=0;
		int total=0;

        int num_matches = (int) matches.size();
		for (int j=0;j < num_matches;j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
            }
			
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate=(double(tp)/actualCorrect);
		newPoint.falseRate=(double(fp)/actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

		dataPoints.push_back(newPoint);
	}

	return dataPoints;
}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	CShape sh = image.Shape();
	Feature f;

	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			double r = image.Pixel(x,y,0);
			double g = image.Pixel(x,y,1);
			double b = image.Pixel(x,y,2);

			if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
				// If the pixel satisfies this meaningless criterion,
				// make it a feature.
				
				f.type = 1;
				f.id += 1;
				f.x = x;
				f.y = y;

				f.data.resize(1);
				f.data[0] = r + g + b;

				features.push_back(f);
			}
		}
	}
}


// Failed attempt to implement HoG, a copy of ComputeHarrisFeatures is below
////void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
////{
////	//Create grayscale image used for Harris detection
////    CFloatImage grayImage = ConvertToGray(image);

////	//Create image to store Harris values
////    CFloatImage harrisImage(image.Shape().width, image.Shape().height, 1);

////	//Create image to store local maximum harris values as 1, other pixels 0
////    CByteImage harrisMaxImage(image.Shape().width, image.Shape().height, 1);


////    // computeHarrisValues puts harris values at each pixel position in harrisImage.
////    computeHarrisValues(grayImage, harrisImage);

////    // Compute max and min value
////    double max_h = -1e100;
////    double min_h = 1e100;
////    for (int y=0; y<harrisImage.Shape().height; y++) {
////        for (int x=0; x<harrisImage.Shape().width; x++) {
////            if (max_h < harrisImage.Pixel(x, y, 0)) {
////                max_h = harrisImage.Pixel(x, y, 0);
////            }
////            if (min_h > harrisImage.Pixel(x, y, 0)) {
////                min_h = harrisImage.Pixel(x, y, 0);
////            }
////        }
////    }

////    // Normalize harrisImage
////    for (int y=0; y<harrisImage.Shape().height; y++) {
////        for (int x=0; x<harrisImage.Shape().width; x++) {
////            harrisImage.Pixel(x, y, 0) = (harrisImage.Pixel(x, y, 0))/max_h;
////        }
////    }

////    // Threshold the harris image and compute local maxima.  You'll need to implement this function.
////    computeLocalMaxima(harrisImage, harrisMaxImage);

////    // Prints out the harris image for debugging purposes
////    CByteImage tmp(harrisImage.Shape());
////    convertToByteImage(harrisImage, tmp);
////    WriteFile(tmp, "harris.tga");


////    // TODO--------------------------------------------------------------------
////    // Loop through feature points in harrisMaxImage and create feature descriptor
////    // for each point above a threshold

////    // Get gradient of harrisImage
////    int w = harrisMaxImage.Shape().width;
////    int h = harrisMaxImage.Shape().height;
////    CFloatImage gradientImageHV(w, h, 3);  // Horizontal and Vertical
////    CFloatImage gradientImageMD(w, h, 3);  // Magnitude and Direction

////    for (int y = 0; y < h-1; y++) {
////        for (int x = 0; x < w-1; x++) {
////            gradientImageHV.Pixel(x, y, 0) = (harrisImage.Pixel(x+1, y, 0)) - (harrisImage.Pixel(x, y, 0)); //h
////            gradientImageHV.Pixel(x, y, 1) = (harrisImage.Pixel(x, y+1, 0)) - (harrisImage.Pixel(x, y, 0)); //v
////            gradientImageMD.Pixel(x, y, 0) = sqrt(pow(gradientImageHV.Pixel(x, y, 0), 2) + pow(gradientImageHV.Pixel(x, y, 1), 2)); //magnitude
////            gradientImageMD.Pixel(x, y, 1) = atan(gradientImageHV.Pixel(x, y, 1) / gradientImageHV.Pixel(x, y, 0)); //direction
////            if (gradientImageMD.Pixel(x, y, 1) < 0) {
////                gradientImageMD.Pixel(x, y, 1) += M_PI/2;
////            }
////        }
////    }

////    int window = 2;
////    int windowWidth = (2*window + 1);
////    int windowSize = pow(windowWidth, 2);

////    double min_angle = 10.0 * M_PI / 180.0;
////    double max_angle = 170.0 * M_PI / 180.0;
////    double step_angle = 20.0 * M_PI / 180.0;
////    int size_hist = (max_angle - min_angle) / step_angle + 1;

////    unsigned int i = 0;
////    for (int y=0; y<harrisMaxImage.Shape().height; y++) {
////        for (int x=0; x<harrisMaxImage.Shape().width; x++) {
////			// Skip over non-maxima
////            if (harrisMaxImage.Pixel(x, y, 0) == 0)
////                continue;

////            // TODO---------------------------------------------------------------------
////		    // Fill in feature with descriptor data here.
////            Feature f;
////            f.type = 1;
////            f.id = i++;
////            f.x = x;
////            f.y = y;

////////////////////////////////////////////////////////////////////////////
///////////// Trying to implement HOG //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////            f.data.resize(size_hist);
////            for (int i = 0; i < size_hist; i++)
////            {
////                f.data[i] = 0;
////            }
////            for (int y_offset = -window; y_offset <= window; y_offset++)
////            {
////                if (y + y_offset < 0 || y + y_offset >= harrisMaxImage.Shape().height) {
////                    continue;
////                }
////                for (int x_offset = -window; x_offset <= window; x_offset++)
////                {
////                    if (x + x_offset < 0 || x + x_offset >= harrisMaxImage.Shape().width) {
////                        continue;
////                    }
////                    double magnitude = gradientPixelMD.Pixel(x + x_offset, y + y_offset, 0);
////                    double angle = gradientPixelMD.Pixel(x + x_offset, y + y_offset, 1);

////                    int i = 0;
////                    for (double angle_h=min_angle; angle_h<=max_angle; angle_h+=step_angle) {
////                        if (angle_h > angle) {
////                            int j = (i-1+size_hist) % (size_hist);
////                            double factor1 = (step_angle - angle_h + angle)/step_angle;
////                            double factor2 = (angle_h - angle)/step_angle;
////                            f.data[i] = f.data[i] + magnitude * factor1;
////                            f.data[j] = f.data[j] + magnitude * factor2;
////                            break;
////                        }
////                        i++;
////                        if (max_angle < angle) {
////                            int j = size_hist - 1;
////                            double factor1 = (-max_angle + angle)/step_angle;
////                            double factor2 = (max_angle + step_angle - angle)/step_angle;
////                            f.data[0] = f.data[0] + magnitude * factor1;
////                            f.data[j] = f.data[j] + magnitude * factor2;
////                            break;
////                        }
////                    }
////                }
////            }
////            max_h = -1e100;
////            min_h = 1e100;
////            double norm = 0;
////            for (int i = 0; i < size_hist; i++)
////            {
////                if (max_h < f.data[i]) {
////                    max_h = f.data[i];
////                }
////                if (min_h > f.data[i]) {
////                    min_h = f.data[i];
////                }
////                norm += pow(f.data[i], 2);
////            }
////            norm = sqrt(norm);

////            for (int i = 0; i < size_hist; i++)
////            {
//////                f.data[i] = (f.data[i] - min_h) / max_h;
////                f.data[i] = (f.data[i]) / norm;
////            }

////            // Add the feature to the list of features
////            features.push_back(f);
////        }
////	}
////}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
    //Create grayscale image used for Harris detection
    CFloatImage grayImage=ConvertToGray(image);

    //Create image to store Harris values
    CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

    //Create image to store local maximum harris values as 1, other pixels 0
    CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);


    //compute Harris values puts harris values at each pixel position in harrisImage.
    //You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage);

    // Compute max and min value
    double max_h = -1e100;
    double min_h = 1e100;
    for (int y=0; y<harrisImage.Shape().height; y++) {
        for (int x=0; x<harrisImage.Shape().width; x++) {
            if (max_h < harrisImage.Pixel(x, y, 0)) {
                max_h = harrisImage.Pixel(x, y, 0);
            }
            if (min_h > harrisImage.Pixel(x, y, 0)) {
                min_h = harrisImage.Pixel(x, y, 0);
            }
        }
    }

    // Normalize harrisImage
    for (int y=0; y<harrisImage.Shape().height; y++) {
        for (int x=0; x<harrisImage.Shape().width; x++) {
            harrisImage.Pixel(x, y, 0) = (harrisImage.Pixel(x, y, 0))/max_h;
        }
    }

    // Threshold the harris image and compute local maxima.  You'll need to implement this function.
    computeLocalMaxima(harrisImage,harrisMaxImage);

    // Prints out the harris image for debugging purposes
    CByteImage tmp(harrisImage.Shape());
    convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");


    // TO DO--------------------------------------------------------------------
    //Loop through feature points in harrisMaxImage and create feature descriptor
    //for each point above a threshold

    int w = harrisMaxImage.Shape().width;
    int h = harrisMaxImage.Shape().height;
    CFloatImage gradientImageMD(w, h, 2);  // Magnitude and Direction

    for (int y = 0; y < h-1; y++) {
        for (int x = 0; x < w-1; x++) {
            double px = harrisImage.Pixel(x+1, y, 0) - harrisImage.Pixel(x, y, 0);;
            double py = harrisImage.Pixel(x, y+1, 0) - harrisImage.Pixel(x, y, 0);;
            gradientImageMD.Pixel(x, y, 1) = atan(py / px); //direction
            if (gradientImageMD.Pixel(x, y, 1) < 0) {
                gradientImageMD.Pixel(x, y, 1) += M_PI/2;
            }
        }
    }


    int window = 3;
    int windowWidth = (2*window + 1);
    int windowSize = pow(windowWidth, 2);

    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
        for (int x=0;x<harrisMaxImage.Shape().width;x++) {

            // Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

            //TO DO---------------------------------------------------------------------
            // Fill in feature with descriptor data here.
            Feature f;
            f.type = 1;
            f.id += 1;
            f.x = x;
            f.y = y;

//            // This is the simplest window descriptor, a 5x5 window around our pixel.
//            // Remember to change window to 2 to see the value used in our examples.
//            // For some reason, sometimes this will still return nan if gradientImageMD loop is not commented
//            f.data.resize(windowSize);
//            for (int y_offset = -window; y_offset <= window; y_offset++)
//            {
//                if (y + y_offset < 0 || y + y_offset >= harrisMaxImage.Shape().height) {
//                    continue;
//                }
//                for (int x_offset = -window; x_offset <= window; x_offset++)
//                {
//                    if (x + x_offset < 0 || x + x_offset >= harrisMaxImage.Shape().width) {
//                        continue;
//                    }
//                    f.data[(x_offset + window) + (windowWidth)*(y_offset + window)] = harrisImage.Pixel(x + x_offset, y + y_offset, 0);
//                }
//            }

            f.data.resize(windowSize*2);
            for (int y_offset = -window; y_offset <= window; y_offset++)
            {
                if (y + y_offset < 0 || y + y_offset >= harrisMaxImage.Shape().height) {
                    continue;
                }
                for (int x_offset = -window; x_offset <= window; x_offset++)
                {
                    if (x + x_offset < 0 || x + x_offset >= harrisMaxImage.Shape().width) {
                        continue;
                    }
                    f.data[(x_offset + window) + (windowWidth)*(y_offset + window)] = harrisImage.Pixel(x + x_offset, y + y_offset, 0);
                    f.data[(x_offset + window) + (windowWidth)*(y_offset + window) + pow(windowWidth, 2)] = gradientImageMD.Pixel(x + x_offset, y + y_offset, 1);
                }
            }

            // Add the feature to the list of features
            features.push_back(f);
        }
    }
}


// TODO---------------------------------------------------------------------
// Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage)
{

    int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;

    CFloatImage gradientImage1(w, h, 1);
    CFloatImage gradientImage2(w, h, 1);
    CFloatImage gradientImage3(w, h, 1);

    for (int y = 0; y < h-1; y++) {
        for (int x = 0; x < w-1; x++) {
            double Ix = (srcImage.Pixel(x+1, y, 0)) - (srcImage.Pixel(x, y, 0));
            double Iy = (srcImage.Pixel(x, y+1, 0)) - (srcImage.Pixel(x, y, 0));
            gradientImage1.Pixel(x, y, 0) = Ix * Ix;
            gradientImage2.Pixel(x, y, 0) = Ix * Iy;
            gradientImage3.Pixel(x, y, 0) = Iy * Iy;
        }
    }

    CFloatImage gaussian(7, 7, 1);
    for (int y = 0; y < 7; y++) {
        for (int x = 0; x < 7; x++) {
            gaussian.Pixel(x, y, 0) = gaussian7x7[x + 7 * y];
        }
    }

    CFloatImage newGradientImage1(w, h, 1);
    CFloatImage newGradientImage2(w, h, 1);
    CFloatImage newGradientImage3(w, h, 1);
    Convolve(gradientImage1, newGradientImage1, gaussian);
    Convolve(gradientImage2, newGradientImage2, gaussian);
    Convolve(gradientImage3, newGradientImage3, gaussian);

//    // Gaussina performed manually
//    for (int y = window; y < h-window; y++) {
//        for (int x = window; x < w-window; x++) {
//            for (int y_offset = -window; y_offset <= window; y_offset++) {
//                int yy = y + y_offset;
//                for (int x_offset = -window; x_offset <= window; x_offset++) {
//                    int xx = x + x_offset;
//                    double w = gaussian[(x_offset+window) + (window*2+1)*(y_offset+window)];
//                    newGradientImage.Pixel(x, y, 0) += gradientImage.Pixel(xx, yy, 0) * w;
//                    newGradientImage.Pixel(x, y, 1) += gradientImage.Pixel(xx, yy, 1) * w;
//                    newGradientImage.Pixel(x, y, 2) += gradientImage.Pixel(xx, yy, 2) * w;
//                }
//            }
//        }
//    }

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            double detH = newGradientImage1.Pixel(x, y, 0)*newGradientImage3.Pixel(x, y, 0) - pow(newGradientImage2.Pixel(x, y, 0), 2);
            double traceH = newGradientImage1.Pixel(x, y, 0) + newGradientImage3.Pixel(x, y, 0);
            if (traceH == 0) {
                traceH = 0.000001;
            }
            harrisImage.Pixel(x, y, 0) = detH / traceH;
        }
    }
}



// TODO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
// You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage, CByteImage &destImage)
{
    double threshold = 0.135;
    int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;
    int window = 1;

    int numPts = 0;
    int numMaxs = 0;
    int numPx = 0;

    for (int y = window; y < h-window; y++)
    {
        for (int x = window; x < w-window; x++)
        {
            numPx++;
            destImage.Pixel(x, y, 0) = 0.0;
            double pixel = srcImage.Pixel(x, y, 0);
            bool localMax = true;
            for (int y_offset = -window; y_offset <= window; y_offset++)
            {
                int yy = y + y_offset;
                if (pixel < srcImage.Pixel(x-1, yy, 0) || pixel < srcImage.Pixel(x, yy, 0) || pixel < srcImage.Pixel(x+1, yy, 0))
                {
                    localMax = false;
                    break;
                }
            }

            if (localMax) {
                numPts++;
                if (pixel > threshold) {
                    destImage.Pixel(x, y, 0) = 1.0;
                    numMaxs += 1;
                }
            }
        }
    }
    printf("\nNum Local Maxima: %d (%.2f%%) (%.3f%%)\n", numMaxs, numMaxs*100.0/numPts, numMaxs*100.0/numPx);
}

// Perform simple feature matching.  This just uses the SSD distance between two feature vectors,
// and matches a feature in the first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

        for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}

        matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest;
		totalScore += matches[i].score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) 
{
    int m = f1.size();
    int n = f2.size();

    matches.resize(m);
    totalScore = 0;

    double d;
    double dBest1;
    int idBest1;
    double dBest2;

    double d_min = 1e100;
    double d_max = -1e100;
    double d_avg = 0;

    for (int i=0; i<m; i++) {
        dBest1 = 1e100;
        idBest1 = 0;
        dBest2 = 1e100;

        for (int j=0; j<n; j++) {
            d = distanceSSD(f1[i].data, f2[j].data);

            if (d < dBest1) {
                dBest2 = dBest1;
                dBest1 = d;
                idBest1 = f2[j].id;
            }
            else if (d < dBest2)
            {
                dBest2 = d;
            }
        }

        matches[i].id1 = f1[i].id;
        matches[i].id2 = idBest1;
        matches[i].score = dBest1/dBest2;
        if (d_min > matches[i].score) {
            d_min = matches[i].score;
        }
        if (d_max < matches[i].score) {
            d_max = matches[i].score;
        }
        d_avg += matches[i].score;
        totalScore += matches[i].score;
    }
    d_avg /= m;

    printf("Min dist: %f\nMax dist: %f\nAvg dist: %f\n", d_min, d_max, d_avg);
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y=0; y<h; y++) {
		for (int x=0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
			}

			index += d;
		}
	}
	
	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<sh.nBands; c++) {
				float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i=0; i<m; i++) {
		dist += pow(v1[i]-v2[i], 2);
	}

	
	return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6]*x + h[7]*y + h[8];

	xNew = (h[0]*x + h[1]*y + h[2]) / d;
	yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc=0;
	double xdiff,ydiff;
	for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff=(results[i].falseRate-results[i-1].falseRate);
		ydiff=(results[i].trueRate-results[i-1].trueRate);
		auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    	    
    }
	return auc;
}

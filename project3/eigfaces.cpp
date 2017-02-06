
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"


EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
    // ----------- TODO #2: compute the coefficients for the face and store in coefficients.
    for (int i = 0; i < getSize(); i++) {
        const Face* u = &((*this)[i]);
        coefficients[i] = 0.0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                coefficients[i] += u->pixel(x, y, 0) * (face.pixel(x, y, 0) - average_face.pixel(x, y, 0));
            }
        }
    }
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
    result.resize(width, height, 1);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            result.pixel(x, y, 0) = average_face.pixel(x, y, 0);
            for (int i = 0; i < getSize(); i++) {
                const Face* u = &((*this)[i]);
                result.pixel(x, y, 0) += coefficients[i] * u->pixel(x, y, 0);
            }
        }
    }
}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.
    Vector coefficients;
    Face result;
    projectFace(face, coefficients);
    constructFace(coefficients, result);
    mse = face.mse(result);
    return mse <= max_reconstructed_mse;
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse.
    Vector coefficients;
    projectFace(face, coefficients);
    mse = coefficients.mse(user_coefficients);
    return mse <= max_coefficients_mse;
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
    Vector coefficients;
    projectFace(face, coefficients);
    for (int i=0; i<users.getSize(); i++) {
        users[i].setMse(coefficients.mse(users[i]));
    }
    users.sort();
}

void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.

    // fill the best_faces vector
    std::list<FacePosition> best_faces;
    for (int i=0; i<n; i++) {
        best_faces.push_back(FacePosition());
    }

    // find the best faces
    for (double scale=min_scale; scale<max_scale+step; scale+=step) {
        Image scaled_img_color((int) img.getWidth()*scale, (int) img.getHeight()*scale, img.getColors());
        img.resample(scaled_img_color);
        Image scaled_img;
        scaled_img_color.grayscale(scaled_img);
        int true_width = (int) width / scale;
        int true_height = (int) height / scale;
        for (int y = 0; y < scaled_img.getHeight() - height + 1; y++) {
            int true_y = (int) y / scale;
            for (int x = 0; x < scaled_img.getWidth() - width + 1; x++) {
                int true_x = (int) x / scale;
                Face tmp_img(width, height);
                double face_mean = 0;
                for (int y_offset=0; y_offset<height; y_offset++) {
                    for (int x_offset=0; x_offset<width; x_offset++) {
                        tmp_img.pixel(x_offset, y_offset, 0) = scaled_img.pixel(x+x_offset, y+y_offset, 0);
                        face_mean += tmp_img.pixel(x_offset, y_offset, 0);
                    }
                }
                face_mean /= (height*width);

                double mse;
                isFace(tmp_img, 1000.0, mse);

                if (mse < best_faces.back().error) {
                    // calculate variance and discard the image if the variance is low
                    double mse_var = 0;
                    for (int y_offset=0; y_offset<height; y_offset++) {
                        for (int x_offset=0; x_offset<width; x_offset++) {
                            mse_var += pow(tmp_img.pixel(x_offset, y_offset, 0) - face_mean, 2);
                        }
                    }
                    if (mse_var < 5e5) {
                        continue;
                    }

                    // see if face overlaps, and if it does, keep only face with best mse
                    int overlap = 0;
                    int pos = 0;
                    double best_pos = 0;
                    double best_mse = DBL_MAX;
                    std::list<int> pos_overlaps;
                    for (std::list<FacePosition>::iterator it=best_faces.begin(); it!=best_faces.end(); it++) {
                        int true_x_it = (int) (*it).x / (*it).scale;
                        int true_y_it = (int) (*it).y / (*it).scale;
                        int true_width_it = (int) width / (*it).scale;
                        int true_height_it = (int) height / (*it).scale;
                        if (true_x_it < true_x+true_width && true_x_it+true_width_it > true_x && true_y_it < true_y+true_height && true_y_it+true_height_it > true_y) {
                            if ((*it).error > mse) {
                                (*it).x = x;
                                (*it).y = y;
                                (*it).scale = scale;
                                (*it).error = mse;
                                (*it).var = mse_var;
                            }
                            if ((*it).error < best_mse) {
                                best_mse = (*it).error;
                                best_pos = pos;
                            }
                            overlap++;
                            pos_overlaps.push_back(pos);
                        }
                        pos++;
                    }
                    if (overlap == 0) {
                        FacePosition best_face;
                        best_face.x = x;
                        best_face.y = y;
                        best_face.scale = scale;
                        best_face.error = mse;
                        best_face.var = mse_var;
                        best_faces.push_back(best_face);
                        best_faces.sort();
                        best_faces.pop_back();
                    } else {
                        // if there were multiple overlaps, remove the overlaps with higher mse
                        if (overlap > 1) {
                            pos = 0;
                            int overlap_pos = pos_overlaps.front();
                            pos_overlaps.pop_front();
                            for (std::list<FacePosition>::iterator it=best_faces.begin(); it!=best_faces.end(); it++) {
                                if (pos == overlap_pos) {
                                    if (pos != best_pos) {
                                        (*it).error = DBL_MAX;
                                    }
                                    if (pos_overlaps.size() == 0) {
                                        break;
                                    }
                                    overlap_pos = pos_overlaps.front();
                                    pos_overlaps.pop_front();
                                }
                                pos++;
                            }
                        }
                        best_faces.sort();
                    }
                }
            }
        }
    }

    if (crop) {
        // return cropped face
        std::list<FacePosition>::iterator it=best_faces.begin();
        int true_x = (int) (*it).x / (*it).scale;
        int true_y = (int) (*it).y / (*it).scale;
        int true_width = (int) width / (*it).scale;
        int true_height = (int) height / (*it).scale;
        img.crop(true_x, true_y, true_x + true_width, true_y + true_height, result);
        std::cout << "Face: scale=" << (*it).scale << " - mse=" << (*it).error << std::endl;
    } else {
        // return original face with squares around the detected faces
        result.resize(img.getWidth(), img.getHeight(), img.getColors());
        img.resample(result);
        int i = 0;
        for (std::list<FacePosition>::iterator it=best_faces.begin(); it!=best_faces.end(); it++) {
            int true_x = (int) (*it).x / (*it).scale;
            int true_y = (int) (*it).y / (*it).scale;
            int true_width = (int) width / (*it).scale;
            int true_height = (int) height / (*it).scale;
            result.line(true_x, true_y, true_x, true_y + true_height, 0, 255, 0);
            result.line(true_x, true_y, true_x + true_width, true_y, 0, 255, 0);
            result.line(true_x + true_width, true_y, true_x + true_width, true_y + true_height, 0, 255, 0);
            result.line(true_x, true_y + true_height, true_x + true_width, true_y + true_height, 0, 255, 0);
            std::cout << "Face" << i << ": scale=" << (*it).scale << " - mse=" << (*it).error << std::endl;
            i++;
        }
    }
}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
    // Negative distances are ok too.

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}




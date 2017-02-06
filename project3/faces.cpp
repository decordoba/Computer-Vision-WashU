
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: faces.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"

#include "jacob.h"
#include <iostream>

Faces::Faces()
:
Array<Face>(),
width(0),
height(0),
vector_size(0)
{
	//empty
}

Faces::Faces(int count, int width, int height)
:
Array<Face>(count),
width(width),
height(height),
vector_size(width*height)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i].resize(width, height, 1);
	}
}

void Faces::load(BinaryFileReader& file)
{
    resize(file.readInt());
	width=file.readInt();
	height=file.readInt();
	vector_size=width*height;
	for (int i=0; i<getSize(); i++) {
		(*this)[i].load(file);
	}
	average_face.load(file);
    std::cout << "Loaded faces from '" << file.getFilename() << "'" << std::endl;

}

void Faces::load(std::string filename)
{
	BinaryFileReader file(filename);
	load(file);
}

void Faces::save(std::string filename) const
{
	BinaryFileWriter file(filename);
	save(file);
}

void Faces::save(BinaryFileWriter& file) const
{
	file.write(getSize());
	file.write(width);
	file.write(height);
	for (int i=0; i<getSize(); i++) {
		(*this)[i].save(file);
	}
	average_face.save(file);
	std::cout << "Saved faces to '" << file.getFilename() << "'" << std::endl;
}

void Faces::output(std::string filepattern) const
{
	for (int i=0; i<getSize(); i++) {
		// normalize for output
		Image out_image;
		(*this)[i].normalize(0.0, 255.0, out_image);
		std::string filename=Functions::filenameNumber(filepattern, i, getSize()-1);
		out_image.saveTarga(filename);
	}
}

void Faces::eigenFaces(EigFaces& results, int n) const
{
	// size the results vector
	results.resize(n);
	results.setHeight(height);
	results.setWidth(width);

    // Change code slightly to allow speed-up and calculate XX' instead of X'X
    bool speed_up = false;
    int matrix_width = vector_size;
    if (speed_up) {
        matrix_width = getSize();
    }

	// allocate matrices
    double **matrix = Jacobi::matrix(1, matrix_width, 1, matrix_width);
    double **eigmatrix = Jacobi::matrix(1, matrix_width, 1, matrix_width);
    double *eigenvec = Jacobi::vector(1, matrix_width);

    // Initialize to zero
    for (int i = 1; i <= matrix_width; i++) {
        for (int j = 0; j <= matrix_width; j++) {
            matrix[i][j] = 0.0;
            eigmatrix[i][j] = 0.0;
        }
        eigenvec[i] = 0.0;
    }

    // --------- TODO #1: fill in your code to prepare a matrix whose eigenvalues and eigenvectors are to be computed.
	// Also be sure you store the average face in results.average_face (A "set" method is provided for this).

    // calculate average face (x_bar)
    Face my_average_face(width, height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int i = 0; i < getSize(); i++) {
                const Face* face = &((*this)[i]);
                my_average_face.pixel(x, y, 0) = my_average_face.pixel(x, y, 0) + face->pixel(x, y, 0);
            }
            my_average_face.pixel(x, y, 0) = my_average_face.pixel(x, y, 0) / getSize();
        }
    }
    results.setAverage(my_average_face);

    // fill matrix (called A in the slides)
    double **faces_diff = Jacobi::matrix(1, getSize(), 1, vector_size);
    for (int i = 0; i < getSize(); i++) {
        const Face* face = &((*this)[i]);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                faces_diff[i + 1][x + y*height + 1] = face->pixel(x, y, 0) - my_average_face.pixel(x, y, 0);
            }
        }
    }
    if (speed_up) {
        for (int i = 1; i <= getSize(); i++) {
            for (int j = 1; j <= getSize(); j++) {
                for (int k = 1; k <= vector_size; k++) {
                    matrix[i][j] += faces_diff[i][k] * faces_diff[j][k];
                }
            }
        }
    } else {
        for (int i = 1; i <= vector_size; i++) {
            for (int j = 1; j <= vector_size; j++) {
                for (int k = 1; k <= getSize(); k++) {
                    matrix[i][j] += faces_diff[k][i] * faces_diff[k][j];
                }
            }
        }
    }

	// find eigenvectors
	int nrot;
    printf("Calculating eigenvectors ...\n");
    Jacobi::jacobi(matrix, matrix_width, eigenvec, eigmatrix, &nrot);
	// sort eigenvectors
	Array<int> ordering;
    for (int i=0; i<ordering.getSize(); i++) {
        std::cout << ordering[i] << std::endl;
    }
    sortEigenvalues(eigenvec, ordering);
    printf("Eigenvectors computed and sorted\n");
	for (int i=0; i<n; i++) {
        for (int k=0; k<matrix_width; k++) {
			results[i][k] = eigmatrix[k+1][ordering[i]+1];
		}
	}
	// free matrices
    Jacobi::free_matrix(matrix, 1, matrix_width, 1, matrix_width);
    Jacobi::free_matrix(eigmatrix, 1, matrix_width, 1, matrix_width);
    Jacobi::free_vector(eigenvec, 1, matrix_width);
}



int Faces::getWidth() const
{
	return width;
}

int Faces::getHeight() const
{
	return height;
}

void Faces::setWidth(int width)
{
	width=width;
	vector_size=width*height;
}

void Faces::setHeight(int height)
{
	height=height;
	vector_size=width*height;
}

void Faces::sortEigenvalues(double *eigenvec, Array<int>& ordering) const
{
	// for now use simple bubble sort
	ordering.resize(vector_size);
	std::list<EigenVectorIndex> list;
	for (int i=0; i<vector_size; i++) {
		EigenVectorIndex e;
		e.eigenvalue=eigenvec[i+1];
		e.index=i;
		list.push_back(e);
	}
	bool change=true;
	list.sort();
	std::list<EigenVectorIndex>::iterator it=list.begin();
	int n=0;
	while (it!=list.end()) {
		ordering[n] = (*it).index;
		it++;
		n++;
	}
}


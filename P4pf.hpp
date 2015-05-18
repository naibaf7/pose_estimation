/*
 * P4pf.h
 *
 *  Created on: Feb 28, 2015
 *      Author: Marco Zorzi
 */

#ifndef P4PF_HPP_
#define P4PF_HPP_

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Eigen>
#include <array>
#include <complex>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/StdVector>

class P4pf {
public:
	P4pf();
	virtual ~P4pf();

	double signOf(double aNumber);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> pointWisePower(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat);

	int getRigidTransform2(Eigen::Matrix<double, 3, Eigen::Dynamic> p1,Eigen::Matrix<double, 3, Eigen::Dynamic> p2,bool bLeftHandSystem, std::vector<Eigen::Matrix<double, 3,3>>* solutions);

	int P4Pf_m(Eigen::Matrix<double, 2, 4> m2D,	Eigen::Matrix<double, 3, 4> M3D, std::vector<double>* focalLengths_m, std::vector<Eigen::Matrix<double, 3,3>>* rotationMatrices, std::vector<Eigen::Matrix<double, 3,1>>* translationVector );

	int P4Pf(Eigen::Matrix<double, 2, 4> m2D,	Eigen::Matrix<double, 3, 4> M3D, std::vector<double>* focalLengths, std::vector<Eigen::Matrix<double, 3,3>>* rotationMatrices, std::vector<Eigen::Matrix<double, 3,1>>* translationVector );

	int fillMatrix(Eigen::Matrix<double,78,88>* matrix, int *arr, double value, int size);

	 int p4pfcode(double glab, double glac, double glad,
			double glbc, double glbd, double glcd,
			double a1, double a2, double b1, double b2,
			double c1, double c2, double d1, double d2, std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>>* solutions);

	 double inline dabs(double a);

	 void GJ(double *A, int rcnt, int ccnt, double tol);

	 void p4pfmex(double *glab, double *a1, double *b1, double *c1, double *d1, Eigen::Matrix<std::complex<double>, 10, 10>  *A);

	 void CalcCoefs(double const *src1, double const *src2, double const *src3, double const *src4, double const *src5, double *dst1);


};

#endif /* P4PF_HPP_ */

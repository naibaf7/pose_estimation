/*
 * Copyright (c) 2011, Laurent Kneip, ETH Zurich
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ETH Zurich nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * P3p.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: Laurent Kneip
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *
 *       Input: featureVectors: 3x3 matrix with UNITARY feature vectors (each column is a vector)
 *              worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
 *              solutions: 3x16 matrix that will contain the solutions
 *                         form: [ 3x1 position(solution1) 3x3 orientation(solution1) 3x1 position(solution2) 3x3 orientation(solution2) ... ]
 *                         the obtained orientation matrices are defined as transforming points from the cam to the world frame
 *      Output: int: 0 if correct execution
 *                  -1 if world points aligned
 */

#include "P3p.hpp"
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>

P3p::P3p() {
}

P3p::~P3p() {
}

int P3p::computePoses(Eigen::Matrix<double, 3, 3> featureVectors,
		Eigen::Matrix<double, 3, 3> worldPoints,
		Eigen::Matrix<double, 3, 16> & solutions) {
	// Extraction of world points

	Eigen::Matrix<double, 3, 1> P1 = worldPoints.col(0);
	Eigen::Matrix<double, 3, 1> P2 = worldPoints.col(1);
	Eigen::Matrix<double, 3, 1> P3 = worldPoints.col(2);

	// Verification that world points are not colinear

	Eigen::Matrix<double, 3, 1> temp1 = P2 - P1;
	Eigen::Matrix<double, 3, 1> temp2 = P3 - P1;

	if ((temp1.cross(temp2)).norm() == 0)
		return -1;

	// Extraction of feature vectors

	Eigen::Matrix<double, 3, 1> f1 = featureVectors.col(0);
	Eigen::Matrix<double, 3, 1> f2 = featureVectors.col(1);
	Eigen::Matrix<double, 3, 1> f3 = featureVectors.col(2);

	// Creation of intermediate camera frame

	Eigen::Matrix<double, 3, 1> e1 = f1;
	Eigen::Matrix<double, 3, 1> e3 = f1.cross(f2);
	e3 = e3 / e3.norm();
	Eigen::Matrix<double, 3, 1> e2 = e3.cross(e1);

	Eigen::Matrix<double, 3, 3> T;
	T.row(0) = e1;
	T.row(1) = e2;
	T.row(2) = e3;

	f3 = T * f3;

	// Reinforce that f3[2] > 0 for having theta in [0;pi]

	if (f3(2) > 0) {
		f1 = featureVectors.col(1);
		f2 = featureVectors.col(0);
		f3 = featureVectors.col(2);

		e1 = f1;
		e3 = f1.cross(f2);
		e3 = e3 / e3.norm();
		e2 = e3.cross(e1);

		T.row(0) = e1;
		T.row(1) = e2;
		T.row(2) = e3;

		f3 = T * f3;

		P1 = worldPoints.col(1);
		P2 = worldPoints.col(0);
		P3 = worldPoints.col(2);
	}

	// Creation of intermediate world frame
	Eigen::Matrix<double, 3, 1> n1 = P2 - P1;
	n1 = n1 / n1.norm();
	Eigen::Matrix<double, 3, 1> n3 = n1.cross(P3 - P1);
	n3 = n3 / n3.norm();
	Eigen::Matrix<double, 3, 1> n2 = n3.cross(n1);

	Eigen::Matrix<double, 3, 3> N;
	N.row(0) = n1;
	N.row(1) = n2;
	N.row(2) = n3;

	// Extraction of known parameters

	P3 = N * (P3 - P1);

	double d_12 = (P2 - P1).norm();
	double f_1 = f3(0) / f3(2);
	double f_2 = f3(1) / f3(2);
	double p_1 = P3(0);
	double p_2 = P3(1);

	double cos_beta = f1.dot(f2);
	double b = 1 / (1 - pow(cos_beta, 2)) - 1;

	if (cos_beta < 0)
		b = -sqrt(b);
	else
		b = sqrt(b);

	// Definition of temporary variables for avoiding multiple computation

	double f_1_pw2 = pow(f_1, 2);
	double f_2_pw2 = pow(f_2, 2);
	double p_1_pw2 = pow(p_1, 2);
	double p_1_pw3 = p_1_pw2 * p_1;
	double p_1_pw4 = p_1_pw3 * p_1;
	double p_2_pw2 = pow(p_2, 2);
	double p_2_pw3 = p_2_pw2 * p_2;
	double p_2_pw4 = p_2_pw3 * p_2;
	double d_12_pw2 = pow(d_12, 2);
	double b_pw2 = pow(b, 2);

	// Computation of factors of 4th degree polynomial

	Eigen::Matrix<double, 5, 1> factors;

	factors(0) = -f_2_pw2 * p_2_pw4 - p_2_pw4 * f_1_pw2 - p_2_pw4;

	factors(1) = 2 * p_2_pw3 * d_12 * b + 2 * f_2_pw2 * p_2_pw3 * d_12 * b
			- 2 * f_2 * p_2_pw3 * f_1 * d_12;

	factors(2) = -f_2_pw2 * p_2_pw2 * p_1_pw2
			- f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2
			- f_2_pw2 * p_2_pw2 * d_12_pw2 + f_2_pw2 * p_2_pw4
			+ p_2_pw4 * f_1_pw2 + 2 * p_1 * p_2_pw2 * d_12
			+ 2 * f_1 * f_2 * p_1 * p_2_pw2 * d_12 * b
			- p_2_pw2 * p_1_pw2 * f_1_pw2 + 2 * p_1 * p_2_pw2 * f_2_pw2 * d_12
			- p_2_pw2 * d_12_pw2 * b_pw2 - 2 * p_1_pw2 * p_2_pw2;

	factors(3) = 2 * p_1_pw2 * p_2 * d_12 * b + 2 * f_2 * p_2_pw3 * f_1 * d_12
			- 2 * f_2_pw2 * p_2_pw3 * d_12 * b - 2 * p_1 * p_2 * d_12_pw2 * b;

	factors(4) = -2 * f_2 * p_2_pw2 * f_1 * p_1 * d_12 * b
			+ f_2_pw2 * p_2_pw2 * d_12_pw2 + 2 * p_1_pw3 * d_12
			- p_1_pw2 * d_12_pw2 + f_2_pw2 * p_2_pw2 * p_1_pw2 - p_1_pw4
			- 2 * f_2_pw2 * p_2_pw2 * p_1 * d_12 + p_2_pw2 * f_1_pw2 * p_1_pw2
			+ f_2_pw2 * p_2_pw2 * d_12_pw2 * b_pw2;

	// Computation of roots

	Eigen::Matrix<double, 4, 1> realRoots;

	this->solveQuartic(factors, realRoots);

	// Backsubstitution of each solution

	for (int i = 0; i < 4; i++) {
		double cot_alpha = (-f_1 * p_1 / f_2 - realRoots(i) * p_2 + d_12 * b)
				/ (-f_1 * realRoots(i) * p_2 / f_2 + p_1 - d_12);

		double cos_theta = realRoots(i);
		double sin_theta = sqrt(1 - pow((double) realRoots(i), 2));
		double sin_alpha = sqrt(1 / (pow(cot_alpha, 2) + 1));
		double cos_alpha = sqrt(1 - pow(sin_alpha, 2));

		if (cot_alpha < 0)
			cos_alpha = -cos_alpha;

		Eigen::Matrix<double, 3, 1> C(
				d_12 * cos_alpha * (sin_alpha * b + cos_alpha),
				cos_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha),
				sin_theta * d_12 * sin_alpha * (sin_alpha * b + cos_alpha));

		C = P1 + N.transpose() * C;

		Eigen::Matrix<double, 3, 3> R;
		R.row(0) << -cos_alpha,
				-sin_alpha * cos_theta, -sin_alpha * sin_theta;
		R.row(1) << sin_alpha,
				-cos_alpha * cos_theta, -cos_alpha * sin_theta;
		R.row(2)  << 0, -sin_theta, cos_theta;

		R = N.transpose() * R.transpose() * T;

		solutions.col(i * 4) = C;
		solutions.col(i * 4 + 1) = R.col(0);
		solutions.col(i * 4 + 2) = R.col(1);
		solutions.col(i * 4 + 3) = R.col(2);
	}

	return 0;
}

int P3p::solveQuartic(Eigen::Matrix<double, 5, 1> factors,
		Eigen::Matrix<double, 4, 1> & realRoots) {
	double A = factors[0];
	double B = factors[1];
	double C = factors[2];
	double D = factors[3];
	double E = factors[4];

	double A_pw2 = A * A;
	double B_pw2 = B * B;
	double A_pw3 = A_pw2 * A;
	double B_pw3 = B_pw2 * B;
	double A_pw4 = A_pw3 * A;
	double B_pw4 = B_pw3 * B;

	double alpha = -3 * B_pw2 / (8 * A_pw2) + C / A;
	double beta = B_pw3 / (8 * A_pw3) - B * C / (2 * A_pw2) + D / A;
	double gamma = -3 * B_pw4 / (256 * A_pw4) + B_pw2 * C / (16 * A_pw3)
			- B * D / (4 * A_pw2) + E / A;

	double alpha_pw2 = alpha * alpha;
	double alpha_pw3 = alpha_pw2 * alpha;

	std::complex<double> P(-alpha_pw2 / 12 - gamma, 0);
	std::complex<double> Q(
			-alpha_pw3 / 108 + alpha * gamma / 3 - pow(beta, 2) / 8, 0);
	std::complex<double> R = -Q / 2.0
			+ sqrt(pow(Q, 2.0) / 4.0 + pow(P, 3.0) / 27.0);

	std::complex<double> U = pow(R, (1.0 / 3.0));
	std::complex<double> y;

	if (U.real() == 0)
		y = -5.0 * alpha / 6.0 - pow(Q, (1.0 / 3.0));
	else
		y = -5.0 * alpha / 6.0 - P / (3.0 * U) + U;

	std::complex<double> w = sqrt(alpha + 2.0 * y);

	std::complex<double> temp;

	temp = -B / (4.0 * A)
			+ 0.5 * (w + sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
	realRoots(0) = temp.real();
	temp = -B / (4.0 * A)
			+ 0.5 * (w - sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
	realRoots(1) = temp.real();
	temp = -B / (4.0 * A)
			+ 0.5 * (-w + sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
	realRoots(2) = temp.real();
	temp = -B / (4.0 * A)
			+ 0.5 * (-w - sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
	realRoots(3) = temp.real();

	return 0;
}

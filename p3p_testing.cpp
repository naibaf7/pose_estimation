/*
 * p3p_testing.cpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Fabian Tschopp
 *
 */
#include <iostream>
#include <chrono>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <chrono>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>
#include "P3p.hpp"
#include "P4pf.hpp"

void testP3p() {

	std::cout << "=================================" << std::endl;
	std::cout << "========== P3p testing ==========" << std::endl;
	std::cout << "=================================" << std::endl;

	P3p p3p;

	Eigen::Matrix<double, 3, 3> points;

	points << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	points << 0.0, 0.0, 1.0, 10.5, 10.5, 10.5, -10.5, -9.5, -9.5;

	Eigen::Matrix<double, 3, 3> features;

	features << -0.049875, -0.049875, 0.049875, 0.049875, -0.049875, -0.049875, 0.99751, 0.99751, 0.99751;

	Eigen::Matrix<double, 3, 16> out;
	Eigen::Matrix<double, 3, 16> golden;

	golden << 1.7739, 0.9917, 0.0064, -0.1281, 0.5665, 1.0000, 0.0003, -0.0066, 0.1272, 0.9993, -0.0019, 0.0372, -0.4751, 0.9952, -0.0049, 0.0977, 0.6782, 0.1283, -0.0356, 0.9911, 0.5227, 0.0066, -0.0059, 1.0000, 0.7136, -0.0366, 0.1435, 0.9890, 0.6218, -0.0977, -0.1019, 0.9900, -9.6385, 0.0018, -0.9994, -0.0361, -9.9416, 0.0003, -1.0000, -0.0059, -11.4235, -0.0072, -0.9896, 0.1434, -8.9839, 0.0051, -0.9948, -0.1019;

	p3p.computePoses(features, points, out);

	std::cout << points << std::endl;
	std::cout << features << std::endl;
	std::cout << out << std::endl;

	double error = (out - golden).norm();

	std::cout << "Error: " << error << std::endl;

	assert(error < 0.15);
}

void testP4pf() {

	std::cout << "=================================" << std::endl;
	std::cout << "========== P4pf testing =========" << std::endl;
	std::cout << "=================================" << std::endl;


	P4pf p4pf;

	Eigen::Matrix<double, 4, 2> feature_vectors;
	Eigen::Matrix<double, 4, 3> world_points;
	std::vector<double> focal_length_solutions;
	std::vector<Eigen::Matrix<double, 3, 3>> rotation_solutions;
	std::vector<Eigen::Matrix<double, 3, 1>> translation_solutions;

	//world_points << -2.5, 5, 2.5, 5, 10, 5, -5, 10, 5, 10, 20, -10;

	world_points << -5, -5, 5, 3, -3, 3, 3, 3, 3, -9, 9, 9;

	feature_vectors << -1, -1, 1, -1, 1, 1, -1, 1;

	p4pf.P4Pf_m(feature_vectors.transpose(), world_points.transpose(), &focal_length_solutions,
			&rotation_solutions, &translation_solutions);

	for (int i = 0; i < focal_length_solutions.size(); ++i) {
		std::cout << focal_length_solutions[i] << std::endl;
		std::cout << translation_solutions[i] << std::endl;
		std::cout << rotation_solutions[i] << std::endl;
	}

}

int main(int argc, char* argv[]) {
	testP3p();
	testP4pf();
}

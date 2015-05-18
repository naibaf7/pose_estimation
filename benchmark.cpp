/*
 * benchmark.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: Fabian Tschopp
 */

#include "benchmark.hpp"
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include "query_loader.hpp"
#include "query_processor.hpp"
#include "import_export.hpp"

namespace pose_estimation {

void Benchmark(QueryLoader &ql, parse_bundler &golden_pb, QueryProcessor *qp) {
	std::ofstream out_file;
	out_file.open("out/benchmark.txt");
	assert(out_file.is_open());

	for (int i = 0; i < ql.GetQueryCount(); ++i) {
		Query golden = LoadGoldenQuery(golden_pb, i,
				"data/Dubrovnik6K/bundle/");

		Query query = ql.GetQuery(i, false);

		std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
		start = std::chrono::high_resolution_clock::now();
		qp->Process(query);
		end = std::chrono::high_resolution_clock::now();

		double elapsed_seconds =
				std::chrono::duration<double>(end - start).count();

		double position_error = (golden.camera_position()
				- query.camera_position()).norm();
		double rotation_error = (golden.camera_rotation()
				- query.camera_rotation()).norm();

		double focal_length_error = (golden.focal_length()
				- query.focal_length());

		std::cout << "================================================="
				<< std::endl;
		std::cout << "Index: " << i + 1 << std::endl;
		std::cout << "Image size: " << query.image_width() << " x "
				<< query.image_height() << std::endl;
		std::cout << "Computation time: " << elapsed_seconds << " s"
				<< std::endl;
		std::cout << "Position error: " << position_error << std::endl;
		std::cout << "Rotation error: " << rotation_error << std::endl;
		std::cout << "Focal length error: " << focal_length_error << std::endl;
		std::cout << "================================================="
				<< std::endl;

		out_file << (i + 1) << "," << std::setprecision(12) << elapsed_seconds
				<< "," << std::setprecision(12) << position_error << ","
				<< std::setprecision(12) << rotation_error << ","
				<< std::setprecision(12) << focal_length_error << std::endl;
		out_file.flush();

	}
	out_file.close();
}

}

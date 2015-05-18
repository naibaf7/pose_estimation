/*
 * pose_estimation.cpp
 *
 *  Created on: Feb 28, 2015
 *      Author: Fabian Tschopp, Marco Zorzi
 */

#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <functional>
#include <iterator>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>

#ifndef __APPLE__
#include <omp.h>
#endif

#include <chrono>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <flann/flann.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include "P3p.hpp"
#include "P4pf.hpp"
#include "SIFT_keypoint.hpp"
#include "parse_bundler.hpp"
#include "bundler_camera.hpp"
#include "query_loader.hpp"
#include "query_processor.hpp"
#include "log.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "import_export.hpp"
#include "benchmark.hpp"

#define DLOG true
#define INFO_DATASET "data/dubrovnik6k.info"
#define GOLDEN_DATASET "data/Dubrovnik6K/bundle/bundle.orig.out"

#define DATABASE_PATH "data/Dubrovnik6K/"
#define QUERY_PATH "data/query/"

using namespace pose_estimation;

int main(int argc, char* argv[]) {

	// For the parallel algorithms, use 4 threads
#ifndef __APPLE__
	omp_set_num_threads(4);
#endif

	// Write some log message
	log("Program started.");

	// Load the INFO file
	parse_bundler pb;

	// Load the golden reference file
	parse_bundler golden_pb;

#pragma omp parallel sections
	{
#pragma omp section
		{
			golden_pb.parse_data(GOLDEN_DATASET, nullptr);
		}
#pragma omp section
		{
			pb.load_from_binary(INFO_DATASET, 1);
		}
	}

	log("Dataset loaded.");

	int processor_type;
	// Select the processor
	std::cout << "0: Basic Processor, 1: Advanced Processor: ";
	std::cin >> processor_type;

	QueryProcessor *qp;

	// Prepare query processor
	if (processor_type == 1) {
		qp = new QueryProcessorAdvanced(pb);
	} else {
		qp = new QueryProcessorBasic(pb);
	}

	// Prepare query loader
	QueryLoader ql(DATABASE_PATH, QUERY_PATH);

	// Load a query file
	int queryindex;
	std::string querypath;
	int mode;
	bool exit = false;

	// User input loop
	while (!exit) {

		std::cout
				<< "0: Exit, 1: Export 3D Model, 2: Process query file, 3: Benchmark: ";
		std::cin >> mode;

		switch (mode) {
		case 0: {
			exit = true;
			break;
		}
		case 1: {
			ExportMesh(golden_pb, "out/", "mesh");
			break;
		}
		case 2: {
			std::cout << "Enter a query index (1 to " << (ql.GetQueryCount())
					<< "): ";

			std::cin >> queryindex;

			int queryidx = std::min(std::max(queryindex - 1, 0),
					ql.GetQueryCount() - 1);

			Query golden = LoadGoldenQuery(golden_pb, queryidx,
					"data/Dubrovnik6K/bundle/");

			Query query = ql.GetQuery(queryidx, true);

			bool success = qp->Process(query);

			std::cout << "Pose Estimation: " << std::endl;
			std::cout << query.focal_length() << std::endl;
			std::cout << query.camera_position() << std::endl;
			std::cout << query.camera_rotation() << std::endl;

			std::cout << "Pose Estimation (Golden): " << std::endl;
			std::cout << golden.focal_length() << std::endl;
			std::cout << golden.camera_position() << std::endl;
			std::cout << golden.camera_rotation() << std::endl;

			ExportCamera(query, "out/", "camera" + std::to_string(queryindex));

			ExportCameraProjectEdges(query, "out/",
					"camera_proj" + std::to_string(queryindex));
			break;
		}
		case 3: {
			Benchmark(ql, golden_pb, qp);
			break;
		}
		default:
			break;
		}
	}

}

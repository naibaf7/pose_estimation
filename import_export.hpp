/*
 * import_export.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: Fabian Tschopp
 */

#ifndef IMPORT_EXPORT_HPP_
#define IMPORT_EXPORT_HPP_

#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include "SIFT_keypoint.hpp"
#include "parse_bundler.hpp"
#include "bundler_camera.hpp"
#include "query_processor.hpp"
#include "query_loader.hpp"

namespace pose_estimation {

void ExportMesh(parse_bundler &parsebundler, std::string path,
		std::string filename);

void ExportCamera(Query &query, std::string path, std::string filename);

Query LoadGoldenQuery(parse_bundler pb, int query_camera_id,
		std::string list_file);

void ExportCameraProjectEdges(Query &query, std::string path,
		std::string filename);

}

#endif /* IMPORT_EXPORT_HPP_ */

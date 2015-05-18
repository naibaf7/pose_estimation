/*
 * import_export.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: Fabian Tschopp
 */

#include "import_export.hpp"
#include <cassert>
#include "query_loader.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include "SIFT_keypoint.hpp"
#include "parse_bundler.hpp"
#include "bundler_camera.hpp"
#include "query_processor.hpp"
#include "pose_utils.hpp"

namespace pose_estimation {

// Export the SfM mesh in PLY format
void ExportMesh(parse_bundler &parsebundler, std::string path,
		std::string filename) {
	std::ofstream out_file;
	out_file.open(path + filename + ".ply");
	assert(out_file.is_open());

	std::vector<feature_3D_info> &sifts = parsebundler.get_feature_infos();

	// Write PLY Header
	out_file << "ply" << std::endl;
	out_file << "format ascii 1.0" << std::endl;
	out_file << "element vertex " << sifts.size() << std::endl;
	out_file << "property float32 x" << std::endl;
	out_file << "property float32 y" << std::endl;
	out_file << "property float32 z" << std::endl;
	out_file << "property uchar red" << std::endl;
	out_file << "property uchar green" << std::endl;
	out_file << "property uchar blue" << std::endl;
	out_file << "end_header" << std::endl;

	out_file << std::fixed;

	for (unsigned int i = 0; i < sifts.size(); ++i) {
		out_file << std::setprecision(16) << sifts[i].point.x << " "
				<< std::setprecision(16) << sifts[i].point.y << " "
				<< std::setprecision(16) << sifts[i].point.z << " "
				<< (int) sifts[i].point.r << " " << (int) sifts[i].point.g
				<< " " << (int) sifts[i].point.b << std::endl;
	}

	out_file.close();
}

// Load the golden model query and pose estimation for comparison
Query LoadGoldenQuery(parse_bundler pb, int query_camera_id,
		std::string list_path) {

	std::ifstream queryfile;
	queryfile.open(list_path + "list.orig.txt");

	std::string delimiter = " ";
	std::string line;

	int counter_query = 0;
	int counter_total = 0;
	int total_camera_id = 0;

	while (getline(queryfile, line)) {

		std::string name = CutStringBefore(line, "/");

		if (name == "query") {
			++counter_query;
		}
		++counter_total;

		if (counter_query == query_camera_id + 1) {
			total_camera_id = counter_total;
			break;
		}
	}

	queryfile.close();

	bundler_camera camera = pb.get_cameras()[total_camera_id - 1];

	// We have camera position and rotation, unlike bundler, so transform:
	Eigen::Matrix<double, 3, 3> xzflip;
	xzflip << -1, 0, 0, 0, 1, 0, 0, 0, -1;
	Eigen::Matrix<double, 3, 3> rotation = camera.rotation;
	Eigen::Matrix<double, 3, 1> position = -rotation.transpose()
			* camera.translation;

	// Convert to our own format
	Query query(nullptr, "", camera.focal_length);
	query.set_camera_position(position);
	rotation = (xzflip * rotation).transpose();
	query.set_camera_rotation(rotation);
	query.set_image_width(camera.width);
	query.set_image_height(camera.height);

	return query;
}

void ExportCameraProjectEdges(Query &query, std::string path,
		std::string filename) {
	std::ofstream out_file_obj;
	out_file_obj.open(path + filename + ".obj");
	assert(out_file_obj.is_open());

	std::vector<Eigen::Matrix<double, 3, 1>> &points_3d = query.fitpoints_3d();

	Eigen::Matrix<double, 3, 1> camera_position = query.camera_position();

	out_file_obj << "o Object" << std::endl;
	out_file_obj << "v " << std::setprecision(16) << camera_position(0) << " "
			<< std::setprecision(16) << camera_position(1) << " "
			<< std::setprecision(16) << camera_position(2) << std::endl;
	out_file_obj << "v " << std::setprecision(16) << camera_position(0) << " "
			<< std::setprecision(16) << camera_position(1) << " "
			<< std::setprecision(16) << camera_position(2) << std::endl;

	for (unsigned int i = 0; i < points_3d.size(); ++i) {
		Eigen::Matrix<double, 3, 1> point_3d = points_3d[i];

		out_file_obj << "v " << std::setprecision(16) << point_3d(0) << " "
				<< std::setprecision(16) << point_3d(1) << " "
				<< std::setprecision(16) << point_3d(2) << std::endl;
	}

	out_file_obj << "s off" << std::endl;
	for (unsigned int i = 0; i < points_3d.size(); ++i) {
		out_file_obj << "f 1 " << i + 3 << " 2" << std::endl;
	}

	out_file_obj.close();
}

// Export the camera pose as OBJ+MTL (wavefront obj with material) and MLP (meshlab project file)
void ExportCamera(Query &query, std::string path, std::string filename) {
	Eigen::Matrix<double, 3, 1> camera_position = query.camera_position();
	Eigen::Matrix<double, 3, 3> camera_rotation = query.camera_rotation();
	Eigen::Matrix<double, 3, 3> xzflip;
	xzflip << -1, 0, 0, 0, 1, 0, 0, 0, -1;
	Eigen::Matrix<double, 3, 3> camera_rotation_xzf = xzflip
			* camera_rotation.transpose();

	std::ofstream out_file_obj;
	std::ofstream out_file_mtl;
	std::ofstream out_file_mlp;
	out_file_obj.open(path + filename + ".obj");
	assert(out_file_obj.is_open());
	out_file_mtl.open(path + filename + ".mtl");
	assert(out_file_mtl.is_open());
	out_file_mlp.open(path + filename + ".mlp");

	// Write obj header
	out_file_obj << "mtllib " << filename << ".mtl" << std::endl;
	out_file_obj << "o Object" << std::endl;

	out_file_obj << std::fixed;

	// Image plane vertices
	for (int i = 0; i < 4; ++i) {
		double x = 0.01 * query.image_width() * std::sin((i * 2 + 1) * M_PI / 4)
				/ sqrt(2);
		double y = 0.01 * query.image_height()
				* std::cos((i * 2 + 1) * M_PI / 4) / sqrt(2);
		double z = query.focal_length() / 100.0;

		// Transform to world coordinates
		Eigen::Matrix<double, 3, 1> xyz;
		xyz << x, y, z;
		xyz = camera_rotation * xyz + camera_position;
		x = xyz(0);
		y = xyz(1);
		z = xyz(2);

		float u = (0.5 + std::sin(-(i * 2 + 1) * M_PI / 4) / sqrt(2.0));
		float v = (0.5 + std::cos((i * 2 + 1) * M_PI / 4) / sqrt(2.0));
		out_file_obj << "v " << std::setprecision(16) << x << " "
				<< std::setprecision(16) << y << " " << std::setprecision(16)
				<< z << std::endl;
		out_file_obj << "vt " << std::setprecision(16) << u << " "
				<< std::setprecision(16) << v << std::endl;
	}

	// Camera box vertices
	out_file_obj << "v " << std::setprecision(16) << camera_position(0) << " "
			<< std::setprecision(16) << camera_position(1) << " "
			<< std::setprecision(16) << camera_position(2) << std::endl;
	for (int i = 0; i < 4; ++i) {
		double x = 0.001 * query.image_width()
				* std::sin((i * 2 + 1) * M_PI / 4) / sqrt(2);
		double y = 0.001 * query.image_height()
				* std::cos((i * 2 + 1) * M_PI / 4) / sqrt(2);
		double z = query.focal_length() / 1000.0;

		// Transform to world coordinates
		Eigen::Matrix<double, 3, 1> xyz;
		xyz << x, y, z;
		xyz = camera_rotation * xyz + camera_position;
		x = xyz(0);
		y = xyz(1);
		z = xyz(2);

		out_file_obj << "v " << std::setprecision(16) << x << " "
				<< std::setprecision(16) << y << " " << std::setprecision(16)
				<< z << std::endl;
	}

	out_file_obj << "g Object_Object_auv" << std::endl;
	out_file_obj << "usemtl Object_auv" << std::endl;
	out_file_obj << "s off" << std::endl;

	// Image plane faces
	out_file_obj << "f 1/1 2/2 3/3" << std::endl;
	out_file_obj << "f 1/1 3/3 4/4" << std::endl;

	// Camera box faces
	out_file_obj << "f 5 6 7" << std::endl;
	out_file_obj << "f 5 7 8" << std::endl;
	out_file_obj << "f 5 8 9" << std::endl;
	out_file_obj << "f 5 6 9" << std::endl;
	out_file_obj << "f 6 7 8" << std::endl;
	out_file_obj << "f 6 8 9" << std::endl;

	// Write mtl
	out_file_mtl << "newmtl Object_auv" << std::endl;
	out_file_mtl << "Ns 100.0" << std::endl;
	out_file_mtl << "d 1.0" << std::endl;
	out_file_mtl << "illum 2" << std::endl;
	out_file_mtl << "Kd 1.0 1.0 1.0" << std::endl;
	out_file_mtl << "Ka 1.0 1.0 1.0" << std::endl;
	out_file_mtl << "Ks 1.0 1.0 1.0" << std::endl;
	out_file_mtl << "Ke 0.0 0.0 0.0" << std::endl;
	out_file_mtl << "map_Kd " << filename << ".jpg" << std::endl;

	// Write mlp

	out_file_mlp << "<!DOCTYPE MeshLabDocument>" << std::endl;
	out_file_mlp << "<MeshLabProject>" << std::endl;
	out_file_mlp << "<RasterGroup>" << std::endl;

	// The camera
	out_file_mlp << "<MLRaster label=" << "\"" << filename + ".jpg" << "\">"
			<< std::endl;
	out_file_mlp << "<VCGCamera ";
	// Translation
	out_file_mlp << "TranslationVector=\"" << std::setprecision(16)
			<< -camera_position(0) << " " << std::setprecision(16)
			<< -camera_position(1) << " " << std::setprecision(16)
			<< -camera_position(2) << "\" ";
	// Various parameters
	out_file_mlp << "LensDistortion=\"0 0\" ";
	out_file_mlp << "ViewportPx=\"" << query.image_width() << " "
			<< query.image_height() << "\" ";
	out_file_mlp << "PixelSizeMm=\"0.001 0.001\" ";
	out_file_mlp << "CenterPx=\"" << query.image_width() / 2 << " "
			<< query.image_height() / 2 << "\" ";
	out_file_mlp << "FocalMm=\"" << query.focal_length() / 1000.0 << "\" ";
	// Rotation matrix
	out_file_mlp << "RotationMatrix=\"" << std::setprecision(16)
			<< camera_rotation_xzf(0, 0) << " " << std::setprecision(16)
			<< camera_rotation_xzf(0, 1) << " " << std::setprecision(16)
			<< camera_rotation_xzf(0, 2) << " " << std::setprecision(16) << 0.0
			<< " " << std::setprecision(16) << camera_rotation_xzf(1, 0) << " "
			<< std::setprecision(16) << camera_rotation_xzf(1, 1) << " "
			<< std::setprecision(16) << camera_rotation_xzf(1, 2) << " "
			<< std::setprecision(16) << 0.0 << " " << std::setprecision(16)
			<< camera_rotation_xzf(2, 0) << " " << std::setprecision(16)
			<< camera_rotation_xzf(2, 1) << " " << std::setprecision(16)
			<< camera_rotation_xzf(2, 2) << " " << std::setprecision(16) << 0.0
			<< " 0.0 0.0 0.0 1.0" << "\" ";
	out_file_mlp << "/>" << std::endl;
	// The image file
	out_file_mlp << "<Plane semantic=\"\" fileName=\"" << filename + ".jpg"
			<< "\"/>" << std::endl;
	out_file_mlp << "</MLRaster>" << std::endl;

	out_file_mlp << "</RasterGroup>" << std::endl;
	out_file_mlp << "</MeshLabProject>" << std::endl;

	out_file_obj.close();
	out_file_mtl.close();
	out_file_mlp.close();

	cv::imwrite(path + filename + ".jpg", query.image());
}

}

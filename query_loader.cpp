/*
 * query_loader.cpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Fabian Tschopp
 */

#include <iostream>
#include <fstream>
#include <string>

#include "parse_bundler.hpp"
#include "bundler_camera.hpp"

#include "query_loader.hpp"
#include "SIFT_keypoint.hpp"
#include "pose_utils.hpp"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace pose_estimation {

QueryLoader::~QueryLoader() {
}

QueryLoader::QueryLoader(std::string database_path, std::string query_path) :
		database_path_(database_path), query_path_(query_path) {

	std::ifstream queryfile;
	queryfile.open(database_path + "list.query.fused.txt");

	std::string delimiter = " ";
	std::string line;

	while (getline(queryfile, line)) {
		double focal_length = 0.0;

		std::vector<std::string> splitline = SplitString(line, " ");

		int width = std::stoi(splitline[1]);
		int height = std::stoi(splitline[2]);

		if (splitline.size() == 4) {
			focal_length = std::stod(splitline[3]);
		}

		std::string name = CutString(splitline[0], "/", ".");

		if (name.length() > 0) {
			Query query(this, name, focal_length);
			query.set_image_width(width);
			query.set_image_height(height);
			queries_.push_back(query);
		}
	}

	queryfile.close();
}

std::string QueryLoader::database_path() {
	return database_path_;
}

std::string QueryLoader::query_path() {
	return query_path_;
}

int QueryLoader::GetQueryCount() {
	return queries_.size();
}

Query& QueryLoader::GetQuery(int index, bool load_image) {
	// Lazy load query image and SIFT when needed
	if(load_image) {
		queries_[index].LoadImage();
	}
	queries_[index].LoadSIFT();
	return queries_[index];
}

Query::~Query() {
}

Query::Query() {
}

Query::Query(QueryLoader *query_loader, std::string query_image_name,
		double focal_length) :
		query_loader_(query_loader), query_image_name_(query_image_name), focal_length_(
				focal_length), image_width_(0), image_height_(0) {
}

bool Query::LoadImage() {

	std::string queryimagepath = query_loader_->query_path() + query_image_name_
			+ ".jpg";

	std::ifstream queryimage(queryimagepath);

	if (queryimage.fail()) {

		std::string img_id = query_image_name_;
		unsigned int img_id_length = img_id.size() + 1;
		while (img_id_length > img_id.size()) {
			img_id_length = img_id.size();
			img_id = CutStringAfter(img_id, "_");
		}

		// Get the webpage containing the flickr image
		std::string wgetcommand =
				"wget -O tmp/tmp.tmp \"http://www.flickr.com/photo_zoom.gne?id="
						+ img_id + "&size=o\"";
		const char* wgetsystemc = wgetcommand.c_str();
		system(wgetsystemc);

		std::ifstream tmpfile;
		tmpfile.open("tmp/tmp.tmp");

		std::string line;

		bool found_image = false;

		while (getline(tmpfile, line)) {
			if (line.find("_o.jpg") != std::string::npos) {
				found_image = true;
				break;
			}
		}

		if (found_image) {
			// Extract the static image path from flickr webpage
			std::string static_image_url = CutString(line, "src=\"", "\">");

			// Get the image itself
			wgetcommand = "wget -O " + queryimagepath + " " + static_image_url;
			wgetsystemc = wgetcommand.c_str();
			system(wgetsystemc);
		} else {
			return false;
		}
	}

	// Actually load the image
	image_ = cv::imread(queryimagepath, CV_LOAD_IMAGE_COLOR);
	//image_width_ = image_.cols;
	//image_height_ = image_.rows;

	return true;
}

bool Query::LoadSIFT() {
	SIFT_loader siftl;

	std::string querypath = query_loader_->query_path() + query_image_name_
			+ ".key";

	std::ifstream querykey(querypath);

	if (querykey.fail()) {
		// Unzip the descriptor if not existing already
		std::string gzipcommand = "gzip -k -d -c "
				+ query_loader_->database_path() + "query/" + query_image_name_
				+ ".key.gz > " + query_loader_->query_path() + query_image_name_
				+ ".key";
		const char* gzipsystemc = gzipcommand.c_str();
		system(gzipsystemc);
	}

	// Load the query descriptors
	const char* querycharpath = querypath.c_str();
	siftl.load_features(querycharpath, LOWE);
	sift_descriptors_ = siftl.get_descriptors();
	sift_keypoints_ = siftl.get_keypoints();
	return true;
}

std::vector<SIFT_keypoint>& Query::sift_keypoints() {
	return sift_keypoints_;
}

std::vector<unsigned char*>& Query::sift_descriptors() {
	return sift_descriptors_;
}

double Query::focal_length() {
	return focal_length_;
}

Eigen::Matrix<double, 3, 1> Query::camera_position() {
	return camera_position_;
}

Eigen::Matrix<double, 3, 3> Query::camera_rotation() {
	return camera_rotation_;
}

void Query::set_fitpoints_2d(
		std::vector<Eigen::Matrix<double, 2, 1>> fitpoints_2d) {
	fitpoints_2d_ = fitpoints_2d;
}
void Query::set_fitpoints_3d(
		std::vector<Eigen::Matrix<double, 3, 1>> fitpoints_3d) {
	fitpoints_3d_ = fitpoints_3d;
}
std::vector<Eigen::Matrix<double, 2, 1>>& Query::fitpoints_2d() {
	return fitpoints_2d_;
}
std::vector<Eigen::Matrix<double, 3, 1>>& Query::fitpoints_3d() {
	return fitpoints_3d_;
}


void Query::set_image_width(int width) {
	image_width_ = width;
}
void Query::set_image_height(int height) {
	image_height_ = height;
}

bool Query::ComputeMatrices() {

	// Camera intrinsic matrix
	camera_matrix_ = Eigen::Matrix<double, 3, 3>::Zero(3, 3);
	camera_matrix_(0, 0) = focal_length_;
	camera_matrix_(1, 1) = focal_length_;
	camera_matrix_(2, 2) = 1.0;
	camera_matrix_(0, 2) = 0.0;//image_width_ / 2;
	camera_matrix_(1, 2) = 0.0;//image_height_ / 2;

	// Inverse matrix for feature preparation
	camera_matrix_inv_ = camera_matrix_.inverse();

	// 4 to 3 conversion identity matrix
	Eigen::Matrix<double, 3, 4> id_mat = Eigen::Matrix<double, 3, 4>::Zero(3,
			4);

	id_mat(0, 0) = 1.0;
	id_mat(1, 1) = 1.0;
	id_mat(2, 2) = 1.0;

	// Camera extrinsic
	Eigen::Matrix<double, 4, 4> rt_mat = Eigen::Matrix<double, 4, 4>::Zero(4,
			4);
	rt_mat.block(0, 0, 3, 3) = camera_rotation_;
	rt_mat.block(0, 3, 3, 1) = camera_position_;
	rt_mat(3, 3) = 1;

	// We have to take the inverse of the RT-matrix because rotation and translation are world->camera.
	proj_matrix_ = camera_matrix_ * id_mat * rt_mat.inverse();

	return true;
}

void Query::set_camera_position(Eigen::Matrix<double, 3, 1> position) {
	camera_position_ = position;
}

void Query::set_camera_rotation(Eigen::Matrix<double, 3, 3> rotation) {
	camera_rotation_ = rotation;
}

void Query::set_focal_length(double focal_length) {
	focal_length_ = focal_length;
}

Eigen::Matrix<double, 3, 3> Query::camera_matrix() {
	return camera_matrix_;
}

Eigen::Matrix<double, 3, 3> Query::camera_matrix_inv() {
	return camera_matrix_inv_;
}

Eigen::Matrix<double, 3, 4> Query::proj_matrix() {
	return proj_matrix_;
}

cv::Mat& Query::image() {
	return image_;
}

int Query::image_width() {
	return image_width_;
}

int Query::image_height() {
	return image_height_;
}

}

/*
 * query_loader.hpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Fabian Tschopp
 */

#ifndef QUERY_LOADER_HPP_
#define QUERY_LOADER_HPP_

#include <string>
#include <vector>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include "SIFT_keypoint.hpp"
#include "parse_bundler.hpp"
#include "bundler_camera.hpp"

namespace pose_estimation {

// Forward declaration
class Query;
class QueryLoader;

// Query class
class Query {
public:
	Query();
	Query(QueryLoader *query_loader, std::string query_image,
			double focal_length);
	virtual ~Query();
	bool LoadImage();
	bool LoadSIFT();
	bool ComputeMatrices();
	double focal_length();
	Eigen::Matrix<double, 3, 1> camera_position();
	Eigen::Matrix<double, 3, 3> camera_rotation();
	void set_camera_position(Eigen::Matrix<double, 3, 1> position);
	void set_camera_rotation(Eigen::Matrix<double, 3, 3> rotation);
	void set_focal_length(double focal_length);
	std::vector<SIFT_keypoint>& sift_keypoints();
	std::vector<unsigned char*>& sift_descriptors();
	Eigen::Matrix<double, 3, 3> camera_matrix();
	Eigen::Matrix<double, 3, 3> camera_matrix_inv();
	Eigen::Matrix<double, 3, 4> proj_matrix();
	cv::Mat& image();
	int image_width();
	int image_height();
	void set_image_width(int width);
	void set_image_height(int height);
	void set_fitpoints_2d(std::vector<Eigen::Matrix<double,2,1>> fitpoints_2d);
	void set_fitpoints_3d(std::vector<Eigen::Matrix<double,3,1>> fitpoints_3d);
	std::vector<Eigen::Matrix<double,2,1>>& fitpoints_2d();
	std::vector<Eigen::Matrix<double,3,1>>& fitpoints_3d();
private:
	QueryLoader *query_loader_;
	std::string query_image_name_;
	double focal_length_;
	Eigen::Matrix<double, 3, 1> camera_position_;
	Eigen::Matrix<double, 3, 3> camera_rotation_;
	Eigen::Matrix<double, 3, 3> camera_matrix_;
	Eigen::Matrix<double, 3, 3> camera_matrix_inv_;
	Eigen::Matrix<double, 3, 4> proj_matrix_;
	std::vector<SIFT_keypoint> sift_keypoints_;
	std::vector<unsigned char*> sift_descriptors_;
	cv::Mat image_;
	int image_width_;
	int image_height_;
	std::vector<Eigen::Matrix<double,2,1>> fitpoints_2d_;
	std::vector<Eigen::Matrix<double,3,1>> fitpoints_3d_;
};

// QueryLoader class
class QueryLoader {
public:
	QueryLoader(std::string database_path, std::string query_path);
	virtual ~QueryLoader();
	std::string database_path();
	std::string query_path();
	Query& GetQuery(int index, bool load_image);
	int GetQueryCount();

private:
	std::string database_path_;
	std::string query_path_;
	std::vector<Query> queries_;
};

}

#endif /* QUERY_LOADER_HPP_ */

/*
 * query_processor_basic.cpp
 *
 *  Created on: Apr 13, 2015
 *      Author: Fabian Tschopp, Marco Zorzi
 */

#include "query_processor.hpp"
#include "query_loader.hpp"
#include "P3p.hpp"
#include "P4pf.hpp"
#include <flann/flann.h>
#include <omp.h>
#include <vector>
#include "pose_utils.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

namespace pose_estimation {

QueryProcessorBasic::QueryProcessorBasic(parse_bundler &parsebundler) :
		QueryProcessor(parsebundler) {

	sifts_ = &(parsebundler_.get_feature_infos());

	features_ = std::vector<float>(sifts_->size() * 128);

	// Collect features and apply averaging
#pragma omp parallel for
	for (unsigned int i = 0; i < sifts_->size(); ++i) {
		for (int j = 0; j < 128; ++j) {
			features_[i * 128 + j] = 0;
		}
		for (unsigned int k = 0; k < ((*sifts_)[i].descriptors.size() / 128);
				++k) {
			for (int j = 0; j < 128; ++j) {
				// 128*k stride, 128 descriptors/stride (j) and i separate features
				features_[i * 128 + j] += (float) (*sifts_)[i].descriptors[j
						+ 128 * k];
			}
		}
		for (int j = 0; j < 128; ++j) {
			features_[i * 128 + j] /= (float) ((*sifts_)[i].descriptors.size()
					/ 128);
		}
	}

	dataset_ = flann::Matrix<float>(&features_[0], sifts_->size(), 128);

	flann_index_ = new flann::Index<flann::L2<float>>(dataset_,
			flann::KDTreeIndexParams(4));
	flann_index_->buildIndex();

}

bool QueryProcessorBasic::Process(Query &query) {

	Query query_copy = query;

	std::vector<unsigned char*> &sifts = query_copy.sift_descriptors();

	std::vector<float> features(sifts.size() * 128);

#pragma omp parallel for
	for (unsigned int i = 0; i < sifts.size(); ++i) {
		for (int j = 0; j < 128; ++j) {
			features[i * 128 + j] = (float) sifts[i][j];
		}
	}

// FLANN kd-tree nearest neighbor search
	int nn = 2;

	flann::Matrix<float> queryset(&features[0], sifts.size(), 128);

	flann::Matrix<int> indices(new int[queryset.rows * nn], queryset.rows, nn);
	flann::Matrix<float> dists(new float[queryset.rows * nn], queryset.rows,
			nn);

	flann_index_->knnSearch(queryset, indices, dists, nn,
			flann::SearchParams(128));

	// Matches that have been correct by SIFT matching
	std::vector<SiftMatch> good_matches;
	// Matches that have been correct by reprojection
	std::vector<SiftMatch> fitted_matches;

// Lowe's ratio test
	float ratio = 0.7;
	for (unsigned int i = 0; i < sifts.size(); ++i) {
		SiftMatch match(i, indices.ptr()[i * nn], dists.ptr()[i * nn]);
		if ((dists.ptr()[i * nn] / dists.ptr()[i * nn + 1] < ratio)) {
			good_matches.push_back(match);
		}
	}

	std::cout << good_matches.size() << std::endl;

	P3p p3p;
	P4pf p4pf;

	bool p4pfmode = query_copy.focal_length() == 0.0;

	// RANSAC
	// Maximum ransac steps
	unsigned int max_ransac_steps = 10000;
	// L2-pixel distance for inliers
	double inlier_eps = 0.5;
	// Fraction of good matches as inliers to stop RANSAC
	unsigned int inlier_divisor = 10;
	// Total amount of good matches as inliers to stop RANSAC
	unsigned int inlier_absolute = 12;
	// Step counter
	unsigned int steps = 0;
	// Hyposet size
	unsigned int hyposet = 0;
	// Hypothesis quality
	double hypoquality = 0;

	// Random sample selector
	std::function<unsigned int()> rand = GetRandomSelector(good_matches.size());

	// Compute initial matrices of the camera (based on image size and focal length)
	if (p4pfmode) {
		// Some "educated guess" on the focal length to get a good inverse camera intrinsic matrix
		query_copy.set_focal_length(query_copy.image_height() / 2);
	}
	query_copy.ComputeMatrices();

// Exit if good model is found or maximum algorithm steps reached
	while (hyposet < good_matches.size() / inlier_divisor
			&& steps < max_ransac_steps && hyposet < inlier_absolute) {

		// Select 3 points for the camera position hypothesis (4 for P4pf)
		std::vector<SiftMatch> selectset;
		for (unsigned int i = 0; i < (p4pfmode == true ? 4 : 3); ++i) {
			bool valid_choice = false;
			while (!valid_choice) {
				SiftMatch select = good_matches[rand()];
				valid_choice = true;
				for (unsigned int j = 0; j < selectset.size(); ++j) {
					SiftMatch opponent = selectset[j];
					if (select.lindex == opponent.lindex) {
						valid_choice = false;
					}
				}
				if (valid_choice) {
					selectset.push_back(select);
				}
			}
		}

		if (p4pfmode) {
			// P4pf case
			Eigen::Matrix<double, 2, 4> feature_vectors;
			Eigen::Matrix<double, 3, 4> world_points;
			std::vector<double> focal_length_solutions;
			std::vector<Eigen::Matrix<double, 3, 3>> rotation_solutions;
			std::vector<Eigen::Matrix<double, 3, 1>> translation_solutions;

			for (int i = 0; i < 4; ++i) {

				SiftMatch match = selectset[i];

				SIFT_keypoint keypoint =
						query_copy.sift_keypoints()[match.lindex];

				point3D point = (*sifts_)[match.mindex].point;

				Eigen::Matrix<double, 3, 1> feature;
				feature
						<< -(keypoint.x
								- (double) (query_copy.image_width()) / 2.0), -(keypoint.y
						- (double) (query_copy.image_height()) / 2.0), query_copy.focal_length();

				feature_vectors.col(i) = feature.block(0, 0, 2, 1);

				// The 3D point
				Eigen::Matrix<double, 3, 1> world_point;
				world_point << point.x, point.y, point.z;
				world_points.col(i) = world_point;
			}

			p4pf.P4Pf_m(feature_vectors, world_points, &focal_length_solutions,
					&rotation_solutions, &translation_solutions);

			// Test the proposed solutions
			for (unsigned int i = 0; i < focal_length_solutions.size(); ++i) {
				double focal_length = focal_length_solutions[i];
				Eigen::Matrix<double, 3, 1> translation =
						translation_solutions[i];
				Eigen::Matrix<double, 3, 3> rotation = rotation_solutions[i];

				Eigen::Matrix<double, 3, 1> camera_position =
						-rotation.transpose() * translation;
				Eigen::Matrix<double, 3, 3> camera_rotation =
						rotation.transpose();

				query_copy.set_focal_length(focal_length);
				query_copy.set_camera_position(camera_position);
				query_copy.set_camera_rotation(camera_rotation);

				// Compute the camera and projection matrix
				query_copy.ComputeMatrices();

				double new_hypoquality = 0;
				std::vector<SiftMatch> fitmatch_vec = TestHypothesis(query_copy,
						good_matches, query_copy.sift_keypoints(), (*sifts_),
						inlier_eps, new_hypoquality);

				unsigned int new_hyposet = fitmatch_vec.size();

				if (new_hypoquality > hypoquality) {
					// New hypothesis was better than the old one
					hypoquality = new_hypoquality;
					hyposet = new_hyposet;
					std::cout << "New hyposet size: " << hyposet
							<< ", with quality " << hypoquality << std::endl;
					// Transfer solution from working copy to original query
					query.set_focal_length(query_copy.focal_length());
					query.set_camera_position(query_copy.camera_position());
					query.set_camera_rotation(query_copy.camera_rotation());
					fitted_matches = fitmatch_vec;
				}
			}

		} else {
			// P3p case
			Eigen::Matrix<double, 3, 3> feature_vectors;
			Eigen::Matrix<double, 3, 3> world_points;
			Eigen::Matrix<double, 3, 16> solutions;

			for (int i = 0; i < 3; ++i) {

				SiftMatch match = selectset[i];

				SIFT_keypoint keypoint =
						query_copy.sift_keypoints()[match.lindex];

				point3D point = (*sifts_)[match.mindex].point;

				Eigen::Matrix<double, 3, 1> feature;
				feature
						<< -(keypoint.x
								- (double) (query_copy.image_width()) / 2.0), -(keypoint.y
						- (double) (query_copy.image_height()) / 2.0), query_copy.focal_length();

				feature.normalize();

				feature_vectors.col(i) = feature;
				Eigen::Matrix<double, 3, 1> world_point;
				world_point << point.x, point.y, point.z;
				world_points.col(i) = world_point;
			}

			p3p.computePoses(feature_vectors, world_points, solutions);

			// Test the 4 proposed solutions
			for (int i = 0; i < 4; ++i) {
				Eigen::Matrix<double, 3, 1> camera_position = solutions.block(0,
						i * 4, 3, 1);
				Eigen::Matrix<double, 3, 3> camera_rotation = solutions.block(0,
						i * 4 + 1, 3, 3);

				query_copy.set_camera_position(camera_position);
				query_copy.set_camera_rotation(camera_rotation);

				// Compute the camera and projection matrix
				query_copy.ComputeMatrices();

				double new_hypoquality = 0;
				std::vector<SiftMatch> fitmatch_vec = TestHypothesis(query_copy,
						good_matches, query_copy.sift_keypoints(), (*sifts_),
						inlier_eps, new_hypoquality);

				unsigned int new_hyposet = fitmatch_vec.size();

				if (new_hypoquality > hypoquality) {
					// New hypothesis was better than the old one
					hypoquality = new_hypoquality;
					hyposet = new_hyposet;
					std::cout << "New hyposet size: " << hyposet
							<< ", with quality " << hypoquality << std::endl;
					// Transfer solution from working copy to original query
					query.set_camera_position(query_copy.camera_position());
					query.set_camera_rotation(query_copy.camera_rotation());
					fitted_matches = fitmatch_vec;
				}
			}
		}
		++steps;
	}

	// Compute final solution on the original query
	query.ComputeMatrices();

	// Refine with Bundle adjustment
	BundleAdjust(query, fitted_matches, query.sift_keypoints(), (*sifts_));

	// Reproject the fitted matches and store in the query
	std::vector<Eigen::Matrix<double, 2, 1>> fitpoints_2d;
	std::vector<Eigen::Matrix<double, 3, 1>> fitpoints_3d;
	Project(query, fitpoints_2d, fitpoints_3d, fitted_matches,
			query.sift_keypoints(), (*sifts_));
	query.set_fitpoints_2d(fitpoints_2d);
	query.set_fitpoints_3d(fitpoints_3d);

	return true;
}

}

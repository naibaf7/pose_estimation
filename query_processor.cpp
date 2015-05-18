/*
 * query_processor.cpp
 *
 *  Created on: Mar 23, 2015
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

SiftMatch::SiftMatch(unsigned int lindex_p, unsigned int mindex_p,
		float distance_p) :
		lindex(lindex_p), mindex(mindex_p), distance(distance_p) {
}

SiftMatch::SiftMatch() :
		lindex(0), mindex(0), distance(0.0) {
}

QueryProcessor::QueryProcessor(parse_bundler &parsebundler) :
		parsebundler_(parsebundler) {

}

bool QueryProcessor::BundleAdjust(Query &query,
		std::vector<SiftMatch> &fitted_matches,
		std::vector<SIFT_keypoint> &sift_keypoints,
		std::vector<feature_3D_info> &sifts) {
	// TODO Levenberg Marquardt for 1 camera pose
	// TODO Write-back directly into query
	// TODO This is optional, if there is enough time to do it.
	return true;
}

bool QueryProcessor::Project(Query &query,
		std::vector<Eigen::Matrix<double, 2, 1>> &fitpoints_2d,
		std::vector<Eigen::Matrix<double, 3, 1>> &fitpoints_3d,
		std::vector<SiftMatch> &matches,
		std::vector<SIFT_keypoint> &sift_keypoints,
		std::vector<feature_3D_info> &sifts) {

	fitpoints_2d.resize(matches.size());
	fitpoints_3d.resize(matches.size());

	Eigen::Matrix<double, 3, 4> proj_matrix = query.proj_matrix();

	// Reproject in parallel
#pragma omp parallel for
	for (unsigned int i = 0; i < matches.size(); ++i) {

		SiftMatch match = matches[i];

		SIFT_keypoint keypoint = sift_keypoints[match.lindex];
		point3D point = sifts[match.mindex].point;

		// Load 3D point
		Eigen::Matrix<double, 4, 1> point3d;
		point3d << point.x, point.y, point.z, 1.0;

		// Load 2D point
		Eigen::Matrix<double, 3, 1> point2d;
		point2d << -(keypoint.x - (double) (query.image_width()) / 2.0), -(keypoint.y
				- (double) (query.image_height()) / 2.0), 1.0;

		// Projection
		Eigen::Matrix<double, 3, 1> point2d_proj = proj_matrix * point3d;

		fitpoints_3d[i] = (point3d.block(0, 0, 3, 1));

		fitpoints_2d[i] = ((point2d_proj / point2d_proj(2)).block(0, 0, 2, 1));

	}

	return true;
}

std::vector<SiftMatch> QueryProcessor::TestHypothesis(Query &query,
		std::vector<SiftMatch> &matches,
		std::vector<SIFT_keypoint> &sift_keypoints,
		std::vector<feature_3D_info> &sifts, double match_eps,
		double &fit_quality) {

	fit_quality = 0;

	std::vector<SiftMatch> fits;

	std::vector<Eigen::Matrix<double, 2, 1>> fitpoints_2d;
	std::vector<Eigen::Matrix<double, 3, 1>> fitpoints_3d;
	Project(query, fitpoints_2d, fitpoints_3d, matches, query.sift_keypoints(),
			(*sifts_));

	for (unsigned int i = 0; i < matches.size(); ++i) {

		SiftMatch match = matches[i];

		SIFT_keypoint keypoint = sift_keypoints[match.lindex];

		// Load 2D point
		Eigen::Matrix<double, 2, 1> point2d;
		point2d << -(keypoint.x - (double) (query.image_width()) / 2.0), -(keypoint.y
				- (double) (query.image_height()) / 2.0);

		// Load 2d reproject point
		Eigen::Matrix<double, 2, 1> point2d_proj;
		point2d_proj = fitpoints_2d[i];

		// Reprojection 2 norm error
		double match_dist = (point2d - point2d_proj).norm();

		if (match_dist < match_eps) {
			fits.push_back(match);
		}
	}

	int width = query.image_width();
	int height = query.image_height();

	Eigen::MatrixXf fit_area(height, width);
	Eigen::MatrixXf match_area(height, width);

	fit_area.setZero();
	match_area.setZero();

	int cover_area = std::ceil((float) (query.image_width()) / 20.0);

#pragma omp parallel for
	for (unsigned int i = 0; i < fits.size(); ++i) {
		SiftMatch &match = fits[i];
		SIFT_keypoint &keypoint = sift_keypoints[match.lindex];
		int x = rint(keypoint.x);
		int y = rint(keypoint.y);

		int start_x = std::max(0, std::min(x - cover_area / 2, width - 1));
		int start_y = std::max(0, std::min(y - cover_area / 2, height - 1));
		int stop_x = std::max(0, std::min(x + cover_area / 2, width - 1));
		int stop_y = std::max(0, std::min(y + cover_area / 2, height - 1));

		if (start_y < stop_y && start_x < stop_x) {
			fit_area.block(start_y, start_x, stop_y - start_y, stop_x - start_x).setOnes();
		}
	}

#pragma omp parallel for
	for (unsigned int i = 0; i < matches.size(); ++i) {
		SiftMatch &match = matches[i];
		SIFT_keypoint &keypoint = sift_keypoints[match.lindex];
		int x = rint(keypoint.x);
		int y = rint(keypoint.y);

		int start_x = std::max(0, std::min(x - cover_area / 2, width - 1));
		int start_y = std::max(0, std::min(y - cover_area / 2, height - 1));
		int stop_x = std::max(0, std::min(x + cover_area / 2, width - 1));
		int stop_y = std::max(0, std::min(y + cover_area / 2, height - 1));

		if (start_y < stop_y && start_x < stop_x) {
			match_area.block(start_y, start_x, stop_y - start_y, stop_x - start_x).setOnes();
		}
	}

	float fit_area_val = fit_area.sum();
	float match_area_val = match_area.sum();

	fit_quality = fit_area_val / match_area_val;

	return fits;
}

}

/*
 * query_processor.hpp
 *
 *  Created on: Mar 23, 2015
 *      Author: Fabian Tschopp, Marco Zorzi
 */

#ifndef QUERY_PROCESSOR_HPP_
#define QUERY_PROCESSOR_HPP_

#include "parse_bundler.hpp"
#include "query_loader.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <flann/flann.h>

namespace pose_estimation {

struct SiftMatch {
	// Index in the match list
	unsigned int lindex;
	// Index of the FLANN matching
	unsigned int mindex;
	// Distance
	float distance;
	// Set of cameras for efficient set intersection
	std::set<unsigned int> camset;
	SiftMatch();
	SiftMatch(unsigned int lindex_p, unsigned int mindex_p, float distance_p);
};

class QueryProcessor {
public:
	QueryProcessor(parse_bundler &parsebundler);
	virtual bool Process(Query &query) = 0;
	std::vector<SiftMatch> TestHypothesis(Query &query,
			std::vector<SiftMatch> &matches,
			std::vector<SIFT_keypoint> &sift_keypoints,
			std::vector<feature_3D_info> &sifts, double match_eps, double &fit_quality);
	bool BundleAdjust(Query &query, std::vector<SiftMatch> &fitted_matches,
			std::vector<SIFT_keypoint> &sift_keypoints,
			std::vector<feature_3D_info> &sifts);
	bool Project(Query &query,
			std::vector<Eigen::Matrix<double, 2, 1>> &fitpoints_2d,
			std::vector<Eigen::Matrix<double, 3, 1>> &fitpoints_3d,
			std::vector<SiftMatch> &matches,
			std::vector<SIFT_keypoint> &sift_keypoints,
			std::vector<feature_3D_info> &sifts);
protected:
	parse_bundler &parsebundler_;
	flann::Index<flann::L2<float>>* flann_index_;
	std::vector<float> features_;
	flann::Matrix<float> dataset_;
	std::vector<feature_3D_info> *sifts_;
};

class QueryProcessorBasic: public QueryProcessor {
public:
	QueryProcessorBasic(parse_bundler &parsebundler);
	bool Process(Query &query);
private:

};

class QueryProcessorAdvanced: public QueryProcessor {
public:
	QueryProcessorAdvanced(parse_bundler &parsebundler);
	bool Process(Query &query);
	bool Backmatching(std::vector<SiftMatch> &good_matches,
			flann::Matrix<float> &queryset);
private:
};

}

#endif /* QUERY_PROCESSOR_HPP_ */

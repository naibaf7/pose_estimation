/*
 * query_processor_advanced.cpp
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
#include <queue>
#include <map>
#include "pose_utils.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>

namespace pose_estimation {

QueryProcessorAdvanced::QueryProcessorAdvanced(parse_bundler &parsebundler) :
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

struct PriorityFeature {
	unsigned int index;
	float priority;
	PriorityFeature(unsigned int index_p, float priority_p);
};

PriorityFeature::PriorityFeature(unsigned int index_p, float priority_p) :
		index(index_p), priority(priority_p) {
}

struct PriorityOrder {
	bool operator()(const PriorityFeature& lhs,
			const PriorityFeature& rhs) const {
		return lhs.priority < rhs.priority;
	}
};

bool QueryProcessorAdvanced::Backmatching(std::vector<SiftMatch> &good_matches,
		flann::Matrix<float> &queryset) {

	// Build a new FLANN KDTree based on image features
	flann::Index<flann::L2<float>> flann_queryindex(queryset,
			flann::KDTreeIndexParams(4));
	flann_queryindex.buildIndex();

	// Number of nearest neighbors to look for
	int nn = 2;
	// Number of backmatches to be achieved
	unsigned int N = 100;
	// Lowe's ratio test value
	float ratio = 0.7;
	// Dynamic priorities booster
	float w = 10;

	flann::Matrix<int> indices(new int[nn], 1, nn);
	flann::Matrix<float> dists(new float[nn], 1, nn);

	std::vector<SiftMatch> backmatches;

	unsigned int t = 0;

	std::priority_queue<PriorityFeature, std::vector<PriorityFeature>,
			PriorityOrder> queue;
	std::priority_queue<PriorityFeature, std::vector<PriorityFeature>,
			PriorityOrder> new_queue;

	// Determine maximum priority from view graph to boost the good matches to the top of the queue
	float max_priority = 0;
	for (unsigned int i = 0; i < sifts_->size(); ++i) {
		int viewlist_size = (*sifts_)[i].view_list.size();
		max_priority = std::max(max_priority, (float) viewlist_size);
	}

	std::cout << "Max priority is " << max_priority << std::endl;

	// Fill in the queue
	for (unsigned int i = 0; i < sifts_->size(); ++i) {
		int viewlist_size = (*sifts_)[i].view_list.size();
		PriorityFeature feature(i, viewlist_size);
		// Boost the good matches to the top of the queue
		for (unsigned int j = 0; j < good_matches.size(); ++j) {
			if (good_matches[j].mindex == i) {
				feature.priority += max_priority;
			}
		}
		queue.push(feature);
	}

	std::map<unsigned int, unsigned int> increase_views;
	//std::vector<unsigned int> increase_views;

	std::cout << "Starting backmatching" << std::endl;
	while (t < 500 * N && backmatches.size() < N) {

		std::vector<float> backfeatures(128);

		// Prioritized picking of a 3D feature point
		PriorityFeature prio_feature = queue.top();
		queue.pop();
#pragma omp parallel for
		for (int i = 0; i < 128; ++i) {
			backfeatures[i] = features_[128 * prio_feature.index + i];
		}

		//std::cout << "FLANN search of feature " << prio_feature.index
		//		<< " with priority " << prio_feature.priority << std::endl;
		flann::Matrix<float> backquery(&backfeatures[0], 1, 128);
		flann_queryindex.knnSearch(backquery, indices, dists, nn,
				flann::SearchParams(128));

		// Ratio test and accept/reject the backmatch
		if ((dists.ptr()[0] / dists.ptr()[1] < ratio)) {
			SiftMatch backmatch(indices.ptr()[0], prio_feature.index,
					dists.ptr()[0]);
			backmatches.push_back(backmatch);
			// Add all views to those that we have to increase in priority
			int viewlist_size = (*sifts_)[prio_feature.index].view_list.size();
			for (int i = 0; i < viewlist_size; ++i) {
				increase_views[(*sifts_)[prio_feature.index].view_list[i].camera]++;
			}
		}

		// Update priority queue according to choices, every 100 steps (excluding the first one)
		if (t != 0 && t % 100 == 0) {
			std::cout << "Backmatching step " << t << ", set size "
					<< backmatches.size() << std::endl;
			// Process every element in the queue (in parallel)
#pragma omp parallel
			{
				while (true) {
					PriorityFeature feature(0, 0);
					bool empty = false;
#pragma omp critical (queuepop)
					{
						if (queue.empty()) {
							empty = true;
						} else {
							feature = queue.top();
							queue.pop();
						}
					}
					if (empty) {
						break;
					}
					int viewlist_size =
							(*sifts_)[feature.index].view_list.size();
					for (std::map<unsigned int, unsigned int>::iterator i =
							increase_views.begin(); i != increase_views.end();
							++i) {
						for (int j = 0; j < viewlist_size; ++j) {
							if ((*sifts_)[feature.index].view_list[j].camera
									== (*i).first) {
								feature.priority += w * (*i).second
										/ (float) viewlist_size;
							}
						}
					}
#pragma omp critical (queuepush)
					{
						new_queue.push(feature);
					}
				}
#pragma omp barrier
			}
			// Employ the updated queue
			std::swap(queue, new_queue);
		}
		++t;
	}
	std::cout << "Backmatching search done, inserting into good matches"
			<< std::endl;

	//good_matches.clear();
	for (unsigned int i = 0; i < backmatches.size(); ++i) {
		// Remove feature duplicates in backmatches
		bool is_duplicate = false;
		for (unsigned int j = 0; j < good_matches.size(); ++j) {
			if (backmatches[i].lindex == good_matches[j].lindex) {
				is_duplicate = true;
			}
		}
		// Copy to good matches
		if (!is_duplicate) {
			// Fill in the camera set
			SiftMatch match = backmatches[i];
			for (unsigned int v = 0;
					v < (*sifts_)[match.mindex].view_list.size(); ++v) {
				match.camset.insert(
						(*sifts_)[match.mindex].view_list[v].camera);
			}
			good_matches.push_back(match);
		}
	}

	std::cout << "Backmatching done, new good matches size is "
			<< good_matches.size() << std::endl;

	// Cleanup and return
	return true;
}

bool QueryProcessorAdvanced::Process(Query &query) {

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
	float ratio = 0.9;
	for (unsigned int i = 0; i < sifts.size(); ++i) {
		if ((dists.ptr()[i * nn] / dists.ptr()[i * nn + 1] < ratio)) {
			// Create a SiftMatch
			SiftMatch match(i, indices.ptr()[i * nn], dists.ptr()[i * nn]);
			// Fill in the camera set
			for (unsigned int v = 0;
					v < (*sifts_)[match.mindex].view_list.size(); ++v) {
				match.camset.insert(
						(*sifts_)[match.mindex].view_list[v].camera);
			}
			good_matches.push_back(match);
		}
	}

	std::cout << " good_matches.size(): " << good_matches.size() << std::endl;

	P3p p3p;
	P4pf p4pf;

	bool p4pfmode = query_copy.focal_length() == 0.0;

	// RANSAC
	// Maximum ransac steps
	unsigned int max_ransac_steps = 100;
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

	// Random probability selector (for improved RANSAC with co-occurence prior)
	std::function<float()> prob = GetRandomProbability();

	// Compute initial matrices of the camera (based on image size and focal length)
	if (p4pfmode) {
		// Some "educated guess" on the focal length to get a good inverse camera intrinsic matrix
		query_copy.set_focal_length(query_copy.image_height() / 2);
	}
	query_copy.ComputeMatrices();

	// Boolean for backmatching
	bool tried_backmatching = false;

	// Exit if good model is found or maximum algorithm steps reached
	while (steps < max_ransac_steps) {

		// Select 3 points for the camera position hypothesis (4 for P4pf)
		//this is what they call K, the subset of of M. M is good_matches
		std::vector<SiftMatch> selectset;
		//Create intersection set

		std::set<unsigned int> intersect;

		// Deterministic upper bound for the selection
		int zero_check_count = 0;
		int zero_check_limit = 200;

		for (unsigned int i = 0; i < (p4pfmode == true ? 4 : 3); ++i) {
			bool valid_choice = false;

			int zero_check = 0;
			while (!valid_choice) {

				std::set<unsigned int> s_camset;
				SiftMatch select;

				// First element of the couple to give to RANSAC
				do {
					// Extracted randomly from good matches
					select = good_matches[rand()];
					//Now we have a valid match
					s_camset = select.camset;
					// We want that more than 5 cameras have seen the first point chosen (i == 0)
					if (i == 0) {
						intersect = s_camset;
					}
				} while (i == 0 && s_camset.size() < 5);

				valid_choice = true;
				for (unsigned int j = 0; j < selectset.size(); ++j) {
					SiftMatch opponent = selectset[j];
					if (select.lindex == opponent.lindex) {
						valid_choice = false;
					}
				}

				// Co-occurence check for i > 0 (points after the first one)
				if (valid_choice && i > 0) {
					//we get here with at least one element in selectset, no duplicates
					std::set<unsigned int> check;
					std::set_intersection(s_camset.begin(), s_camset.end(),
							intersect.begin(), intersect.end(),
							std::inserter(check, check.begin()));

#ifdef DBG
					std::cout << "s_camset.size()= " << s_camset.size() << std::endl;
					std::cout << "intersect.size()= " << intersect.size() << std::endl;
					std::cout << "check.size()= " << check.size() << std::endl;
#endif

					if (check.size() == 0
							&& zero_check_count <= zero_check_limit) {
						zero_check++;
						valid_choice = false;
					} else {
						float ratio = (float) (check.size())
								/ (std::min(intersect.size(), s_camset.size()));
						// 75% acceptance probability pre-multiplier at check.size() == k (==5)
						float k = 5;
						ratio *= 1.0 / (1.0 + exp(-(check.size() / k)));
						float temp = prob();

						// Accept everything to get to a deterministic upper bound solution
						if (zero_check_count > zero_check_limit) {
							ratio = 1.0;
						}

#ifdef DBG
						std::cout << "GetRandomProbability()= " << temp << std::endl;
						std::cout << "ratio= " << ratio << std::endl;
#endif

						if (ratio > temp) {
#ifdef DBG
							std::cout << "Point accepted!" << std::endl;
#endif

							valid_choice = true;
							intersect = check;
						} else {
							valid_choice = false;
#ifdef DBG
							std::cout << "Point not accepted" << std::endl;
#endif
						}
					}
				}

				if (valid_choice) {
					selectset.push_back(select);
				}

				if (zero_check > 30) {
#ifdef DBG
					std::cout << "wooops, Zerocheck failed" << std::endl;
#endif
					// Reset loop, new first point gets chosen
					i = -1;
					selectset.clear();
					intersect.clear();
					valid_choice = true;
					zero_check_count++;
				}
			}
#ifdef DBG
			std::cout << "\nEnded While valid choice\n\n";
#endif
		}
#ifdef DBG
		std::cout << "Selectset.size(): " << selectset.size() << std::endl;
#endif

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

		// Increase the RANSAC step counter
		++steps;

		// If the RANSAC steps are used up and no solution is found, we can try backmatching and reset the RANSAC steps
		if (!tried_backmatching
				&& hyposet < good_matches.size() / inlier_divisor
				&& hyposet < inlier_absolute && steps == max_ransac_steps) {

			// Augment the good matches by doing backmatching:
			Backmatching(good_matches, queryset);

			// Replace the random function, because good_matches size has changed
			rand = GetRandomSelector(good_matches.size());

			// Reset RANSAC steps
			steps = 0;

			// Reset hyposet
			// hyposet = 0;
			// hypoquality = 0;

			// Backmatching applied, mark as true
			tried_backmatching = true;
		}
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


/*
 * pose_utils.cpp
 *
 *  Created on: Mar 26, 2015
 *      Author: Fabian Tschopp
 */

#include <sys/time.h>
#include "pose_utils.hpp"
#include <random>
#include <functional>

namespace pose_estimation {

std::vector<std::string> SplitString(std::string input, std::string delimiter) {
	std::vector<std::string> output;

	std::string cutline = input;

	size_t pos = 0;
	std::string token;
	while ((pos = cutline.find(delimiter)) != std::string::npos) {
		token = cutline.substr(0, pos);
		cutline.erase(0, pos + delimiter.length());
		output.push_back(token);
	}

	output.push_back(cutline);

	return output;
}

std::string CutString(std::string input, std::string delim_from,
		std::string delim_to) {
	size_t from = input.find(delim_from);
	size_t to = input.find(delim_to);
	return input.substr(from + delim_from.length(),
			to - from - delim_from.length());
}

std::string CutStringAfter(std::string input, std::string delim_from) {
	size_t from = input.find(delim_from);
	if (from == std::string::npos) {
		return input;
	}
	return input.substr(from + delim_from.length(),
			input.length() - from - delim_from.length());
}

std::string CutStringBefore(std::string input, std::string delim_to) {
	size_t to = input.find(delim_to);
	if (to == std::string::npos) {
		return input;
	}
	return input.substr(0, to);
}

std::function<unsigned int()> GetRandomSelector(unsigned int set_size) {
	struct timeval start_time;
	gettimeofday(&start_time, NULL);
#ifdef __APPLE__
	std::seed_seq seq { (int)(start_time.tv_sec), (int)(start_time.tv_usec) };
#else
	std::seed_seq seq { start_time.tv_sec, start_time.tv_usec };
#endif
	std::mt19937_64 generator(seq);
	std::uniform_int_distribution<unsigned int> distribution(0, set_size - 1);
	std::function<unsigned int()> selector = std::bind(distribution, generator);
	return selector;
}

std::function<float()> GetRandomProbability() {
	struct timeval start_time;
	gettimeofday(&start_time, NULL);
#ifdef __APPLE__
	std::seed_seq seq { (int)(start_time.tv_sec), (int)(start_time.tv_usec) };
#else
	std::seed_seq seq { start_time.tv_sec, start_time.tv_usec };
#endif
	std::mt19937_64 generator(seq);
	std::uniform_real_distribution<float> distribution(0, 1);
	std::function<float()> probability = std::bind(distribution, generator);
	return probability;
}

}

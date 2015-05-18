/*
 * pose_utils.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Fabian Tschopp
 */

#ifndef POSE_UTILS_HPP_
#define POSE_UTILS_HPP_

#include <string>
#include <vector>
#include <functional>

namespace pose_estimation {

std::vector<std::string> SplitString(std::string input, std::string delimiter);

std::string CutString(std::string input, std::string delim_from,
		std::string delim_to);

std::string CutStringAfter(std::string input, std::string delim_from);

std::string CutStringBefore(std::string input, std::string delim_to);

std::function<unsigned int()> GetRandomSelector(unsigned int set_size);

std::function<float()> GetRandomProbability();

}

#endif /* POSE_UTILS_HPP_ */

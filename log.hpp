/*
 * log.hpp
 *
 *  Created on: Mar 21, 2015
 *      Author: Fabian Tschopp
 */

#ifndef LOG_HPP_
#define LOG_HPP_

#include <iostream>


#ifndef DLOG
#define DLOG true
#endif


void log(std::string info) {
	if (DLOG) {
		struct timeval tv;
		gettimeofday(&tv, 0);
		time_t now = tv.tv_sec;
		struct tm * now_tm = localtime(&now);
		char buffer[40];
		strftime(buffer, sizeof(buffer), "%Y/%m/%d %H:%M:%S", now_tm);
		std::cout << "[" << buffer << "] -- " << info << std::endl;
	}
}



#endif /* LOG_HPP_ */

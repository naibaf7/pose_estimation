/*
 * benchmark.hpp
 *
 *  Created on: Mar 31, 2015
 *      Author: Fabian Tschopp
 */

#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include "benchmark.hpp"
#include "query_loader.hpp"
#include "query_processor.hpp"

namespace pose_estimation {

void Benchmark(QueryLoader &ql, parse_bundler &golden_pb, QueryProcessor *qp);

}

#endif /* BENCHMARK_HPP_ */

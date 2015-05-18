/*===========================================================================*\
 *                                                                           *
 *                            ACG Localizer                                  *
 *      Copyright (C) 2011-2012 by Computer Graphics Group, RWTH Aachen      *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of ACG Localizer                                       *
 *                                                                           *
 *  ACG Localizer is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  ACG Localizer is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with ACG Localizer.  If not, see <http://www.gnu.org/licenses/>.   *
 *                                                                           *
 \*===========================================================================*/

#include "bundler_camera.hpp"

bundler_camera::bundler_camera() {
	rotation.setIdentity();
	focal_length = 1.0;
	kappa_1 = kappa_2 = 0.0;
	translation[0] = translation[1] = translation[2] = 0.0;
	width = height = 0;
}

bundler_camera::~bundler_camera() {
}

//-----------------------------------------------------------------------------

Eigen::Matrix<float, 2, 1> bundler_camera::project_f(
		const Eigen::Matrix<double, 3, 1> &p) const {
	Eigen::Matrix<double, 3, 1> tmp;
	tmp = rotation * p;
	tmp += translation;
	tmp[2] *= -1.0;
	for (int i = 0; i < 2; ++i)
		tmp[i] /= tmp[2];

	double length_sqr = pow((double) tmp[0], 2.0) + pow((double) tmp[1], 2.0);
	// compute radial distortion term
	double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;
	return Eigen::Matrix<float, 2, 1>(focal_length * r * tmp[0],
			focal_length * r * tmp[1]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<float, 2, 1> bundler_camera::project_f_undist(
		const Eigen::Matrix<double, 3, 1> &p) const {
	Eigen::Matrix<double, 3, 1> tmp;
	tmp = rotation * p;
	tmp += translation;
	tmp[2] *= -1.0;
	for (int i = 0; i < 2; ++i)
		tmp[i] /= tmp[2];

	return Eigen::Matrix<float, 2, 1>(focal_length * tmp[0],
			focal_length * tmp[1]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<double, 2, 1> bundler_camera::project_d(
		const Eigen::Matrix<double, 3, 1> &p) const {
	Eigen::Matrix<double, 3, 1> tmp;
	tmp = rotation * p;
	tmp += translation;
	tmp[2] *= -1.0;
	for (int i = 0; i < 2; ++i)
		tmp[i] /= tmp[2];

	double length_sqr = pow((double) tmp[0], 2.0) + pow((double) tmp[1], 2.0);
	// compute radial distortion term
	double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;
	return Eigen::Matrix<double, 2, 1>(focal_length * r * tmp[0],
			focal_length * r * tmp[1]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<double, 2, 1> bundler_camera::project_d_undist(
		const Eigen::Matrix<double, 3, 1> &p) const {
	Eigen::Matrix<double, 3, 1> tmp;
	tmp = rotation * p;
	tmp += translation;
	tmp[2] *= -1.0;
	for (int i = 0; i < 2; ++i)
		tmp[i] /= tmp[2];

	return Eigen::Matrix<double, 2, 1>(focal_length * tmp[0],
			focal_length * tmp[1]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<float, 3, 1> bundler_camera::get_cam_global_vec_f(
		const Eigen::Matrix<double, 3, 1> &_local) const {
	Eigen::Matrix<double, 3, 1> g_pos;
	g_pos = _local.transpose() * rotation;
	return Eigen::Matrix<float, 3, 1>(g_pos[0], g_pos[1], g_pos[2]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<double, 3, 1> bundler_camera::get_cam_global_vec_d(
		const Eigen::Matrix<double, 3, 1> &_local) const {
	Eigen::Matrix<double, 3, 1> g_pos;
	g_pos = _local.transpose() * rotation;
	return g_pos;
}

//-----------------------------------------------------------------------------

Eigen::Matrix<float, 3, 1> bundler_camera::get_cam_position_f() const {
	Eigen::Matrix<double, 3, 1> g_pos;
	g_pos = translation.transpose() * rotation;
	return Eigen::Matrix<float, 3, 1>(-g_pos[0], -g_pos[1], -g_pos[2]);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<double, 3, 1> bundler_camera::get_cam_position_d() const {
	Eigen::Matrix<double, 3, 1> g_pos;
	g_pos = translation.transpose() * rotation;
	return Eigen::Matrix<double, 3, 1>(-g_pos[0], -g_pos[1], -g_pos[2]);
}

//-----------------------------------------------------------------------------

double bundler_camera::compute_reprojection_error_f(
		Eigen::Matrix<float, 3, 1> &p, float x, float y) const {
	double error = 0.0;

	Eigen::Matrix<double, 3, 1> tmp, tmp2(p[0], p[1], p[2]);
	tmp = tmp2.transpose() * rotation;
	tmp += translation;
	tmp[2] *= -1.0;
	for (int i = 0; i < 2; ++i)
		tmp[i] /= tmp[2];

	double length_sqr = pow((double) tmp[0], 2.0) + pow((double) tmp[1], 2.0);
	// compute radial distortion term
	double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;

	Eigen::Matrix<double, 2, 1> proj_point(focal_length * r * (double) tmp[0],
			focal_length * r * (double) tmp[1]);

	error += (x - proj_point[0]) * (x - proj_point[0]);
	error += (y - proj_point[1]) * (y - proj_point[1]);

//   std::cout << proj_point << " vs ( " << x << " , " << y << " ) " << std::endl;

	return sqrt(error);
}

//-----------------------------------------------------------------------------

Eigen::Matrix<double, 3, 1> bundler_camera::to_cam_coords_d(
		const Eigen::Matrix<double, 3, 1> &p) const {
	Eigen::Matrix<double, 3, 1> tmp;
	tmp = rotation * p;
	tmp += translation;

	return tmp;
}


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

#ifndef BUNDLER_CAMERA_HH
#define BUNDLER_CAMERA_HH

/**
 * Class to model a Camera in a Bundler reconstruction.
 *
 * author: Torsten Sattler (tsattler@cs.rwth-aachen.de)
 * date  : 09-26-2011
 **/

#include <stdint.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>

/**
 * Class for storing a camera consisting of a rotation matrix,
 * a translation vector, a focal length value and two parameters 
 * accounting for radial distortion
 *
 * Author: Torsten Sattler <tsattler@cs.rwth-aachen.de>
 * Date: 08-02-2010 
 **/

class bundler_camera {
public:
	//! constructor
	bundler_camera();

	//! destructor
	~bundler_camera();

	//! project a 3D point to 2D float coordinates (centered around the principal point of the camera)
	Eigen::Matrix<float, 2, 1> project_f(
			const Eigen::Matrix<double, 3, 1> &p) const;

	//! project a 3D point to 2D float coordinates (centered around the principal point of the camera), without distorting it
	Eigen::Matrix<float, 2, 1> project_f_undist(
			const Eigen::Matrix<double, 3, 1> &p) const;

	//! project a 3D point to 2D double coordinates (centered around the principal point of the camera)
	Eigen::Matrix<double, 2, 1> project_d(
			const Eigen::Matrix<double, 3, 1> &p) const;

	//! project a 3D point to 2D double coordinates (centered around the principal point of the camera), without distorting it
	Eigen::Matrix<double, 2, 1> project_d_undist(
			const Eigen::Matrix<double, 3, 1> &p) const;

	//! compute the 3D direction of the camera in world space
	Eigen::Matrix<float, 3, 1> get_cam_global_vec_f(
			const Eigen::Matrix<double, 3, 1> &_local) const;

	//! compute the 3D direction of the camera in world space
	Eigen::Matrix<double, 3, 1> get_cam_global_vec_d(
			const Eigen::Matrix<double, 3, 1> &_local) const;

	//! compute the 3D position of the camera in world space
	Eigen::Matrix<float, 3, 1> get_cam_position_f() const;

	//! compute the 3D position of the camera in world space
	Eigen::Matrix<double, 3, 1> get_cam_position_d() const;

	//! compute the reprojection error of a 3D point and a 2D point
	double compute_reprojection_error_f(Eigen::Matrix<float, 3, 1> &p, float x,
			float y) const;

	//! compute the coordinates of the point in camera coordinates (remember, the camera looks down the -z axis!)
	Eigen::Matrix<double, 3, 1> to_cam_coords_d(
			const Eigen::Matrix<double, 3, 1> &p) const;

	Eigen::Matrix<double, 3, 3> rotation;
	Eigen::Matrix<double, 3, 1> translation;
	double focal_length;
	double kappa_1, kappa_2;
	uint32_t id;
	int width, height;

};

#endif 

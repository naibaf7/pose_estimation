/*
* p3p_testing.cpp
*
*  Created on: Mar 21, 2015
*      Author: Marco Zorzi
*
*/
#include <iostream>
#include <chrono>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <array>
#include <string>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <chrono>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>
#include "P4pf.hpp"

#define DBG_TEST 0
#define OUT 1

void testP4p_m() {
	
	/*********** BEGIN main funciton ************/
	/*
	% P4P + unknown focal length testing app
	% given a set of 4x 2D<->3D correspondences, calculate camera pose and the
	% camera focal length.
	%
	% by Martin Bujnak, (c)apr2008
	%
	
	% ground truth focal, image resolution [-0.5;0.5]x[-0.5;0.5]
	*/
	
	double fgt = 1.5;
	
	// general scene
	Eigen::Matrix<double,3,4> M;
	M << -3.33639834336120, -23.35549638285873, -13.18941519576778,  6.43164913914748,
	0.65948286155096,  -4.90376715918747,   1.17103701629876,  0.14580433383203,
	-8.46658219501120,  -3.99876939947909, -3.02248927651177,  -22.16086539862748;
	
	Eigen::Matrix<double,2,4> m;
	m << 0.11009888473695,   0.39776592879400,   0.28752996253253,   -0.05017617729940,
	0.03882059658299,  -0.17303453640632,  -0.05791310109713,    0.19297848817239;
	
	// ground truth orientation + position
	Eigen::Matrix<double,3,3> Rgt;
	Rgt << -0.94382954756954,   0.33043272406750,                  0,
	0.27314119331857,   0.78018522420862,  -0.56276540799790,
	-0.18595610677571,  -0.53115462041845,  -0.82661665574857;
	
	Eigen::Matrix<double,3,1> Tgt;
	Tgt << 2.85817409844554,
	-2.17296255562889,
	77.54246130780075;
	
	std::cout << M << std::endl;
	std::cout << m << std::endl;
	std::cout << Rgt << std::endl;
	std::cout << Tgt << std::endl;

	P4pf object;
	
	std::vector<double> focalLengthsResults;
	std::vector<Eigen::Matrix<double, 3,3>> rotationMatrices;
	std::vector<Eigen::Matrix<double, 3,1>> translationResult;
	
	if (DBG_TEST) std::cout << "\n\n focalLengthsResults.size()= " << focalLengthsResults.size() <<"\n";
	if (DBG_TEST) std::cout << "\n\n rotationMatrices.size()= " << rotationMatrices.size() <<"\n";
	if (DBG_TEST) std::cout << "\n\n translationResult.size()= " << translationResult.size() <<"\n";
	
	
	
	
	if (OUT) std::cout <<     "******************************************** LEAVING TESTING *********************************************\n\n";
	std::chrono::milliseconds ms_start = std::chrono::duration_cast< std::chrono::milliseconds >(
		std::chrono::high_resolution_clock::now().time_since_epoch()	);
	
	
	int checkError = object.P4Pf_m(m, M, &focalLengthsResults, &rotationMatrices, &translationResult);
	
	std::chrono::milliseconds ms_stop = std::chrono::duration_cast< std::chrono::milliseconds >(
		std::chrono::high_resolution_clock::now().time_since_epoch()	);
	
	if (OUT)  std::cout << "***************************************** RETURNED TO TESTING ******************************************\n\n";
	if (OUT)  std::cout << "Completed in "<< (ms_stop.count()-ms_start.count()) <<" millisecond(ms)"  << "with checkError = "<< checkError <<"\n\n";
	
	if (DBG_TEST) std::cout << "\n focalLengthsResults.size()= " << focalLengthsResults.size() <<"\n";
	if (DBG_TEST) std::cout << "\n rotationMatrices.size()= " << rotationMatrices.size() <<"\n";
	if (DBG_TEST) std::cout << "\n translationResult.size()= " << translationResult.size() <<"\n";
	
	
	if (DBG_TEST) std::cout << "\n focalLengthsResults:\n";
	for (int i = 0; i< focalLengthsResults.size(); i++)
	{
		if (DBG_TEST) std::cout << focalLengthsResults.at(i) <<"\n";
	}
	if (DBG_TEST) std::cout<<"\n";
	
	if (DBG_TEST) std::cout << "\n rotationMatrices:\n";
	for (int i = 0; i< rotationMatrices.size(); i++)
	{
		if (DBG_TEST) std::cout << "R("<<i<<"):\n"<< rotationMatrices.at(i) <<"\n";
	}
	if (DBG_TEST) std::cout<<"\n";
	
	if (DBG_TEST) std::cout << "\n translationResult:\n";
	for (int i = 0; i< translationResult.size(); i++)
	{
		if (DBG_TEST) std::cout << "T("<<i<<"):\n"<< translationResult.at(i) <<"\n";
	}
	if (DBG_TEST)std::cout<<"\n";
	
	
	Eigen::Matrix<double,3,3> Rrel = Eigen::Matrix<double,3,3>::Zero(3,3);
	Eigen::Matrix<double,3,1> Trel = Eigen::Matrix<double,3,1>::Zero(3,1);
	
	double dangle = -1;
	for (int i = 0; i< focalLengthsResults.size(); i++)
	{
		
		if (DBG_TEST) std::cout << "\n i = " << i<<"\n\n";
		
		if (DBG_TEST) std::cout << "\n Rgt:\n" << Rgt <<"\n\n";
		if (DBG_TEST) std::cout << "\n Rgt.inverse():\n" << Rgt.inverse() <<"\n\n";
		if (DBG_TEST) std::cout << "\n rotationMatrices.at(i):\n" << rotationMatrices.at(i) <<"\n\n";
		
		
		Rrel = Rgt.inverse()*rotationMatrices.at(i);
		if (DBG_TEST) std::cout << "\n Rrel:\n" << Rrel <<"\n\n";
		
		Trel = Tgt - translationResult.at(i);
		if (DBG_TEST) std::cout << "\n Trel:\n" << Trel <<"\n\n";
		
		if (DBG_TEST) std::cout << "Rrel.trace():\n" << Rrel.trace() <<"\n\n";
		if (DBG_TEST) std::cout << "Rrel.trace() -1:\n" << Rrel.trace() -1 <<"\n\n";
		if (DBG_TEST) std::cout << "(Rrel.trace() -1)/2:\n" << (Rrel.trace() -1)/2 <<"\n\n";
		if (DBG_TEST) std::cout << "acos(1):\n" << acos(1) <<"\n\n";
		if (DBG_TEST) std::cout << "acos(1-0.000000000001):\n" << acos(1-0.000000000001) <<"\n\n";
		if (DBG_TEST) std::cout << "acos( (Rrel.trace() -1)/2 ):\n" << std::acos( (Rrel.trace() -1)/2 ) <<"\n\n";
		if (DBG_TEST) std::cout << "std::norm(acos( (Rrel.trace() -1)/2 )):\n" << std::norm(std::acos( (Rrel.trace() -1)/2 )) <<"\n\n";
		if (DBG_TEST) std::cout << "sqrt(std::norm(acos( (Rrel.trace() -1)/2 ))):\n" << std::sqrt(std::norm(std::acos( (Rrel.trace() -1)/2 ))) <<"\n\n";
		if (DBG_TEST) std::cout << "(sqrt(std::norm(acos( (Rrel.trace() -1)/2 ))) )/M_PI*180:\n" << (std::sqrt(std::norm(std::acos( (Rrel.trace() -1)/2 ))) )/M_PI*180 <<"\n\n";
		
		dangle = (
			sqrt(
				std::norm(
					acos(
						(Rrel.trace() -1)/2
						)
					)
				)
			)/M_PI*180;
		if (DBG_TEST) std::cout << "dangle:\n" << dangle <<"\n\n";
		if (DBG_TEST) std::cout << "focalLengthsResults.at(i)=    " << focalLengthsResults.at(i) <<"\n\n";
		
		if (OUT) std::cout << "focal err= " << fgt-focalLengthsResults.at(i) << "\nRotation err= " << dangle << "\nTranslation err= "<<Trel.norm() << "\n\n";
	}
	
	
}

void testP4p() {
	
	Eigen::Matrix<double,3,4> M;
	M << -3.33639834336120, -23.35549638285873, -13.18941519576778,  6.43164913914748,
	0.65948286155096,  -4.90376715918747,   1.17103701629876,  0.14580433383203,
	-8.46658219501120,  -3.99876939947909, -3.02248927651177,  -22.16086539862748;
	
	Eigen::Matrix<double,2,4> m;
	m << 0.11009888473695,   0.39776592879400,   0.28752996253253,   -0.05017617729940,
	0.03882059658299,  -0.17303453640632,  -0.05791310109713,    0.19297848817239;
	
	// ground truth orientation + position
	Eigen::Matrix<double,3,3> Rgt;
	Rgt << -0.94382954756954,   0.33043272406750,                  0,
	0.27314119331857,   0.78018522420862,  -0.56276540799790,
	-0.18595610677571,  -0.53115462041845,  -0.82661665574857;
	
	Eigen::Matrix<double,3,1> Tgt;
	Tgt << 2.85817409844554,
	-2.17296255562889,
	77.54246130780075;
	
	P4pf object;
	
	double fgt = 1.5;
	
	
	std::vector<double> focalLengthsResults2;
	
	std::vector<Eigen::Matrix<double, 3,3>> rotationMatrices2;
	
	std::vector<Eigen::Matrix<double, 3,1>> translationResult2;
	
	if (DBG_TEST) std::cout << "\n\n focalLengthsResults2.size()= " << focalLengthsResults2.size() <<"\n";
	if (DBG_TEST) std::cout << "\n\n rotationMatrices2.size()= " << rotationMatrices2.size() <<"\n";
	if (DBG_TEST) std::cout << "\n\n translationResult2.size()= " << translationResult2.size() <<"\n";
	
	
	
	if (OUT) std::cout <<     "******************************************** LEAVING TESTING *********************************************\n\n";
	
	std::chrono::milliseconds ms_start2 = std::chrono::duration_cast< std::chrono::milliseconds >(
		std::chrono::high_resolution_clock::now().time_since_epoch()	);
	
	
	int checkError2 = object.P4Pf(m, M, &focalLengthsResults2, &rotationMatrices2, &translationResult2);
	
	std::chrono::milliseconds ms_stop2 = std::chrono::duration_cast< std::chrono::milliseconds >(
		std::chrono::high_resolution_clock::now().time_since_epoch()	);
	
	if (OUT)  std::cout << "***************************************** RETURNED TO TESTING ******************************************\n\n";
	
	if (OUT)  std::cout << "Completed in "<< (ms_stop2.count()-ms_start2.count()) <<" millisecond(ms)"  << "with checkError2 = "<< checkError2 <<"\n\n";
	
	if (DBG_TEST) std::cout << "\n focalLengthsResults2.size()= " << focalLengthsResults2.size() <<"\n";
	if (DBG_TEST) std::cout << "\n rotationMatrices2.size()= " << rotationMatrices2.size() <<"\n";
	if (DBG_TEST) std::cout << "\n translationResult2.size()= " << translationResult2.size() <<"\n";
	
	
	if (DBG_TEST) std::cout << "\n focalLengthsResults2:\n";
	for (int i = 0; i< focalLengthsResults2.size(); i++)
	{
		if (DBG_TEST) std::cout << focalLengthsResults2.at(i) <<"\n";
	}
	if (DBG_TEST) std::cout<<"\n";
	
	if (DBG_TEST) std::cout << "\n rotationMatrices:\n";
	for (int i = 0; i< rotationMatrices2.size(); i++)
	{
		if (DBG_TEST) std::cout << "R("<<i<<"):\n"<< rotationMatrices2.at(i) <<"\n";
	}
	if (DBG_TEST) std::cout<<"\n";
	
	if (DBG_TEST) std::cout << "\n translationResult:\n";
	for (int i = 0; i< translationResult2.size(); i++)
	{
		if (DBG_TEST) std::cout << "T("<<i<<"):\n"<< translationResult2.at(i) <<"\n";
	}
	if (DBG_TEST) std::cout<<"\n";
	
	
	Eigen::Matrix<double,3,3> Rrel2 = Eigen::Matrix<double,3,3>::Zero(3,3);
	Eigen::Matrix<double,3,1> Trel2 = Eigen::Matrix<double,3,1>::Zero(3,1);
	
	double dangle2 = -1;
	for (int i = 0; i< focalLengthsResults2.size(); i++)
	{
		
		if (DBG_TEST) std::cout << "\n i = " << i<<"\n\n";
		
		if (DBG_TEST) std::cout << "\n Rgt:\n" << Rgt <<"\n\n";
		if (DBG_TEST) std::cout << "\n Rgt.inverse():\n" << Rgt.inverse() <<"\n\n";
		if (DBG_TEST) std::cout << "\n rotationMatrices.at(i):\n" << rotationMatrices2.at(i) <<"\n\n";
		
		
		Rrel2 = Rgt.inverse()*rotationMatrices2.at(i);
		if (DBG_TEST) std::cout << "\n Rrel2:\n" << Rrel2 <<"\n\n";
		
		Trel2 = Tgt - translationResult2.at(i);
		if (DBG_TEST) std::cout << "\n Trel2:\n" << Trel2 <<"\n\n";
		
		dangle2 = (
			sqrt(
				std::norm(
					acos(
						(Rrel2.trace() -1)/2
						)
					)
				)
			)/M_PI*180;
		
		if (DBG_TEST) std::cout << "dangle:\n" << dangle2 <<"\n\n";
		if (DBG_TEST) std::cout << "focalLengthsResults.at(i)=    " << focalLengthsResults2.at(i) <<"\n\n";
		if (OUT) std::cout << "focal err= " << fgt-focalLengthsResults2.at(i) << "\nRotation err= " << dangle2 << "\nTranslation err= "<<Trel2.norm() << "\n\n";
	}
}

int main(int argc, char* argv[]) {
	if (OUT) std::cout << "\n*********************************************  P4PF_M BEGIN  *********************************************\n\n";
	testP4p_m();
	if (OUT) std::cout << "*********************************************  P4PF_M END  ***********************************************\n\n";
	if (OUT) std::cout << "*********************************************  P4PF BEGIN  ***********************************************\n\n";
	testP4p();
	if (OUT) std::cout << "*********************************************  P4PF END  *************************************************\n\n";
	if (OUT) std::cout << "************************************** P4P runned successfully *******************************************\n\n";
	return 0;
}

/*
* P4pf.cpp
*
*  Created on: Mar 15, 2015
*      Author: Marco Zorzi
*/

#include "P4pf.hpp"			// this includes the header. All other includes are included there.



#define DBG 0				// debug flag for P4pf function
#define DBG_M 0				// debug flag for P4pf_m function
#define DBG_PWP 0			// debug flag for pointWisePower function
#define DBG_FILL 0			// debug flag for fillMatrix function
#define DBG_P4PFCODE 0		// debug flag for p4pfcode function
#define DBG_TR 0			// debug flag for getRigidTransform2 function


P4pf::P4pf() {
}

P4pf::~P4pf() {
}

/* INLINE FUNCTION FOR CALCULATING THE SIGN OF A NUMBER
this was made in order to manipulate matrices, since it needs to 
return ONLY 1, -1 or 0.
INPUT: a double number
OUTPUT: 1  if input number is >0
-1 if input number is <0
0  if input number is =0
*/
double inline P4pf::signOf(double aNumber)
{
	if (aNumber > 0)
		return 1;
	else if (aNumber <0)
		return -1;
	else
		return 0;
}

/* FUNCTION THAT COMPUTES THE 2ND POINTWISE POWER OF A MATRIX OF ANY SIZE 
This was not found in the standard library for Dynamic Matrices. 

INPUT: a Eigen Matrix of generic size
OUTPUT: the given matrix with its element squared
*/
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P4pf::pointWisePower(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat)
{																																																												{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result = (Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(mat.rows(),mat.cols()));
	
	if (DBG_PWP) std::cout << "mat,rows()" << mat.rows() << "mat cols() " << mat.cols() <<" mat size() " <<mat.size() <<"\n";
	
	for(int i=0; i<mat.rows(); i++)
	{
		if (DBG_PWP) std::cout << "i= " << i <<"\n";
		for(int j=0; j<mat.cols(); j++)
		{
			if (DBG_PWP) std::cout << "j= " << j <<"\n";
			double temp = mat(i,j);
			if (DBG_PWP) std::cout << "temp= " << temp <<"\n";
			result(i,j)= temp*temp;
			if (DBG_PWP) std::cout << "temp*temp= " << temp*temp <<"\n";
		}
	}
	return result;
}																																																											}
/*********************************************/
/*********** BEGIN first funciton ************/
/*********************************************/

/* COMMENTS ON THE GIVEN MATLAB CODE
% find rigid transformation (rotation translation) from p1->p2 given
% 3 or more points in the 3D space
%
% based on : K.Arun, T.Huangs, D.Blostein. Least-squares
%            fitting of two 3D point sets IEEE PAMI 1987
%
% by Martin Bujnak, nov2007
%
%function [R t] = GetRigidTransform(p1, p2, bLeftHandSystem)
% p1 - 3xN matrix with reference 3D points
% p2 - 3xN matrix with target 3D points
% bLeftHandSystem - camera coordinate system
*/

/* 
INPUT: 	- p1 and p2, two matrices of 3XN size.
- boolean for the camera coordiante system
- a vector of matrices for the output
OUTPUT:	- Outputs a vector of Eigen 3x3 matrices with these elements: 1)the first element is the rotation matrix
2) the second element is the translation vector
- the int return is used as error check: if 0 correct, if -1 an error was found

NOTE: In order to increase efficiency it was preferred to fit a 3X1 vector in a 3X3 matrix in order to allow Eigen library
to pre-allocate the memory. This way, the array doesn't need to be of Eigen 3XDynamic size. Nevertheless,
the second element of the vector is just a 3x1 vector and the other elements are zeros.
*/
int P4pf::getRigidTransform2(Eigen::Matrix<double, 3, Eigen::Dynamic> p1,Eigen::Matrix<double, 3, Eigen::Dynamic> p2,bool bLeftHandSystem, std::vector<Eigen::Matrix<double, 3,3>>* solutions)
{
	
	//prints out the element received
	if (DBG_TR) std::cout << "p1\n" <<p1 << "\n\n";
	if (DBG_TR) std::cout << "p2\n"<< p2 << "\n\n";
	
	//creation of means matrices
	Eigen::Matrix<double, 3, 1> p1mean; 
	Eigen::Matrix<double, 3, 1> p2mean ;
	//Means calculations
	p1mean = p1.rowwise().sum();
	p1mean = p1mean/ 4;
	p2mean = p2.rowwise().sum();
	p2mean = p2mean/ 4;
	// Prints the means after the calculations. Used to verify rowwise sum of Eigen Library
	if (DBG_TR) std::cout << "p1mean: \n" << p1mean << "\n\n";
	if (DBG_TR) std::cout << "p2mean:\n" << p2mean << "\n\n";
	
	//Replicates the mean in 4 columns and subtracts from the received matrix
	p1 = p1 - p1mean.replicate(1, 4);
	p2 = p2 - p2mean.replicate(1, 4);
	// Print updated matrices to verify replicate function of Eigen Library
	if (DBG_TR) std::cout << "p1\n" <<p1 << "\n\n";
	if (DBG_TR) std::cout << "p2\n"<< p2 << "\n\n";
	
	//normalize to unit size
	// Used very spacy coding for better understanding of the operations order.
	// Saved in a temporary matrix to reduce complexity.
	Eigen::Matrix<double, 1, Eigen::Dynamic> temp1 =
	(
		(
			(
				(
					(
						this->pointWisePower(p1)	// matrix gets squared pointwise
						)
					).colwise()									
					).sum()												// colwise sum
					).cwiseSqrt()												// pointwise square root
					).cwiseInverse();													// Pointwise inversion
					
					// final operations and then save in u1  		
					Eigen::Matrix<double, 3, Eigen::Dynamic> u1 = p1.cwiseProduct(temp1.replicate(3,1)); //replicate horizontal vector 3 times
					//same operations for u2
					Eigen::Matrix<double, 1, Eigen::Dynamic> temp2 =
					(
						(
							(
								(
									(
										this->pointWisePower(p2)
										)
									).colwise()
									).sum()
									).cwiseSqrt()
									).cwiseInverse();
									Eigen::Matrix<double, 3, Eigen::Dynamic> u2 = p2.cwiseProduct(temp2.replicate(3,1));
									
									// Prints u1 and u2 after the calculations
									if (DBG_TR) std::cout << "u1\n" <<u1 << "\n\n";
									if (DBG_TR) std::cout << "u2\n"<< u2 << "\n\n";
									
									// calcolations for rotation
									Eigen::Matrix<double, 3, 3> C = u2 * u1.transpose();
									// prints C for verifications
									if (DBG_TR) std::cout << "C created. Here it is\n" << C << "\n";
									
									// calculations of singular value decomposition using EigenSolver JacboiSVD
									// Std outputs used to verify the progression of the process
									if (DBG_TR) std::cout << "Ready for SVD\n";
									
									// The two flags "ComputeFullU and ComputeFullV tells the solver to calculate the full squared matrices.
									Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
									if (DBG_TR) std::cout << "Made SVD\n";
									
									// Jacobi SVD doesn't have flags for giving the singular values as a matrix.
									// this needs to be done
									Eigen::Matrix<double, 3, 1> S_diag = svd.singularValues();
									Eigen::Matrix<double, 3, 3> S = Eigen::Matrix<double, 3, 3>::Zero(3,3);
									
									// diagonal matrix creation
									if (DBG_TR) std::cout << "Got S\n";
									for (int i = 0; i < 3; i++)
									{
										S(i,i) = S_diag(i);
									}
									// Printed S to verify if it's correct
									if (DBG_TR) std::cout << "S\n" << S << "\n\n";
									
									// Calculated U and printed to verify if it's correct
									Eigen::Matrix<double, 3, 3> U = svd.matrixU();
									if (DBG_TR) std::cout << "U\n" << U << "\n\n";
									
									// Calculated V and printed to verify if it's correct
									Eigen::Matrix<double, 3, 3> V = svd.matrixV();
									if (DBG_TR) std::cout << "V\n" << V << "\n\n";
									
									// Cast to double is required here because the access to a Matrix element somehow
									// returns a matrix subset.
									S(0,0) = (double) this->signOf((double) S(0,0));
									S(1,1) = (double) this->signOf((double) S(1,1));
									
									if (bLeftHandSystem)
										S(2,2) = -(double) this->signOf((double) (U*V.transpose()).determinant());
									else
										S(2,2) = (double) this->signOf((double) (U*V.transpose()).determinant());
									
									//Check of S
									if (DBG_TR) std::cout << "filled S\n";
									if (DBG_TR) std::cout << "S\n" << S << "\n\n";
									
									// Final calculations of R and T, than verified by printing them
									Eigen::Matrix<double, 3, 3> R = U*S*(V.transpose());
									Eigen::Matrix<double, 3, 1> t = -R*p1mean + p2mean;
									if (DBG_TR) std::cout << "R and t\n";
									if (DBG_TR) std::cout << "R\n" << R << "\n\n";
									if (DBG_TR) std::cout << "t\n" << t << "\n\n";
									
									// as per function introduction, here the translation vector (3X1) gets stuffed
									// in the first column of a 3x3 matrix and then pushed back on the solution vector
									Eigen::Matrix<double, 3, 3> t_out = Eigen::Matrix<double, 3, 3>::Zero(3,3);
									for(int i = 0; i<3;i++)
									{
										t_out(i,0)=t(i);
									}
									//Verification of the output
									if (DBG_TR) std::cout << "T out ready\n";
									if (DBG_TR) std::cout << "t_out\n" << t_out << "\n\n";
									
									//solution given back
									solutions->push_back(R);
									solutions->push_back(t_out);
									return 0;
}

/*********************************************/
/************ END first funciton *************/
/*********************************************/

/*********************************************/
/*********** BEGIN second funciton ***********/
/*********************************************/

/* Comments on the original Matlab code
% P4P + unknown focal length
% given a set of 4x 2D<->3D correspondences, calculate camera pose and the
% camera focal length.
%
% by Martin Bujnak, (c)apr2008
%
%
% Please refer to the following paper, when using this code :
%
%      Bujnak, M., Kukelova, Z., and Pajdla, T. A general solution to the p4p
%      problem for camera with unknown focal length. CVPR 2008, Anchorage,
%      Alaska, USA, June 2008
%
%
% function [f R t] = P4Pf_m(m2D, M3D)
%
% input:
%
%  m2D - 2x4 matrix with 4x2D measuremets
%  M3D - 3x4 matrix with corresponding 3D points
%
% output:
%
%  f - vector with N focal lengths
%  R - 3x3xN matrix with N rotation matrices
%  t - 3xN matrix with N translation matrices
%
%  following equation holds for each solution
%
%      lambda*m2D = diag([f(i) f(i) 1])*[R(:,:,i) t(:,i)] * M3D
*/
/* 
INPUT: 	- m2D, 2x4 matrix
- M3D, 3x4 matrix
- a vector of double for returning focal lenghts
- a vector of matrices 3x3 for the rotation matrices
- a vector of Eigen vectors 3x1 for the translation vectors
OUTPUT:	- The given vectors filled
- the int return is used as error check: if 0 correct, if -1 an error was found
*/ 
int P4pf::P4Pf_m(Eigen::Matrix<double, 2, 4> m2D,	Eigen::Matrix<double, 3, 4> M3D, std::vector<double>* focalLengths_mResults, std::vector<Eigen::Matrix<double, 3,3>>* rotationMatrices, std::vector<Eigen::Matrix<double, 3,1>>* translationResult )
{
	//check data arrived
	if (DBG_M) std::cout << "\n**************** ENTERING P4PF_M ****************\n";
	if (DBG_M) std::cout << "\nThis is M3D\n" << M3D << "\n";
	if (DBG_M) std::cout << "\nThis is m2D\n" << m2D << "\n";
	
	// shift 3D data so that variance = sqrt(2), mean = 0
	Eigen::Matrix<double, 3, 1> mean3d = M3D.rowwise().sum() /4;
	M3D = M3D - mean3d.replicate(1, 4);
	
	// Created this matrix for calculation of var. M3d needs to stay as it is. 
	Eigen::Matrix<double, 3, 4> M3Dbis = Eigen::Matrix<double, 3, 4>::Zero(3,4);
	
	// Check if all the calculations give appropriate results.
	// Several bugs solved here
	if (DBG_M) std::cout << "\nThis is M3D\n" << M3D << "\n";
	if (DBG_M) std::cout << "\nThis is M3D^2\n" << M3Dbis << "\n";
	if (DBG_M) std::cout << "\nThis is M3D colwise \n" << M3Dbis.colwise().sum() << "\n";
	if (DBG_M) std::cout << "\nThis is M3D colwise + sqrt \n" << M3Dbis.colwise().sum().cwiseSqrt() << "\n";
	if (DBG_M) std::cout << "\nThis is M3D colwise + sqrt + rowise sum\n" << M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum() << "\n";
	if (DBG_M) std::cout << "\nThis is M3D colwise + sqrt + rowise sum at 0\n" << M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum()(0) << "\n";
	
	// Calculation of var in several step to reduce errors.
	M3Dbis = this->pointWisePower(M3D);
	double var = M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum()(0);
	var = var/4;

	if (DBG_M) std::cout << "\nThis is var = " << var << "\n";

	M3D = (1/var)*M3D;
	
	// Check of var and of M3D modified
	if (DBG_M) std::cout << "\nThis is var = " << var << "\n";
	if (DBG_M) std::cout << "\nThis is M3D final, after it was elaborated with var\n" << M3D << "\n";
	
	// scale 2D data
	double var2d = this->pointWisePower(m2D).colwise().sum().cwiseSqrt().rowwise().sum()(0);
	var2d = var2d/4;
	if (DBG_M) std::cout << "\n var2d extracted from m2D = " << var2d;


	m2D = (1/var2d)*m2D;
	
	// Checking of the results
	if (DBG_M) std::cout << "\n var2d extracted from m2D = " << var2d <<"\n";
	if (DBG_M) std::cout << "\nThis is m2D final, after it was elaborated with var2D\n" << m2D << "\n";
	
	//coefficients of 5 reduced polynomials in these monomials mon
	double glab = this->pointWisePower( (M3D.col(0) - M3D.col(1)) ).colwise().sum()(0);
	double glac = this->pointWisePower( (M3D.col(0) - M3D.col(2)) ).colwise().sum()(0);
	double glad = this->pointWisePower( (M3D.col(0) - M3D.col(3)) ).colwise().sum()(0);
	double glbc = this->pointWisePower( (M3D.col(1) - M3D.col(2)) ).colwise().sum()(0);
	double glbd = this->pointWisePower( (M3D.col(1) - M3D.col(3)) ).colwise().sum()(0);
	double glcd = this->pointWisePower( (M3D.col(2) - M3D.col(3)) ).colwise().sum()(0);
	
	// Consistecy check of the coefficients
	if (glbc*glbd*glcd*glab*glac*glad < 1e-15)
	{
		return -1;
	}
	if (DBG_M) std::cout << "\n Passed the coefficient check, the are not zero\n";
	
	// Creating array for solutions and calling the function
	std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> solutions;
	if (DBG_M) std::cout << "solutions.size()= " << solutions.size() << "\n";
	
	if (DBG_M) std::cout << "\n**************** LEAVING P4PF_M ****************\n";
	int errorCheck= this->P4pf::p4pfcode(glab, glac, glad, glbc, glbd, glcd, (double)m2D(0,0),(double) m2D(1,0),
		(double)m2D(0,1),(double) m2D(1,1),(double) m2D(0,2), (double)m2D(1,2),
		(double) m2D(0,3),(double) m2D(1,3), &solutions);
	if (DBG_M) std::cout << "\n**************** RE-ENTERED P4PF_M after P4PFCODE ****************\n";
	
	//Checking if there are errors
	if (DBG_M) std::cout << "\n called p4pfcode. errorCheck = " << errorCheck <<".\n";
	if (errorCheck == -1)
	{
		return -1;
	}
	// Check of solution consistency
	if (DBG_M) std::cout << "solutions.size()= " << solutions.size() << "\n";
	
	// Unpacking of the solutions
	// This was tricky because in matlab they are created in an opposite way.
	Eigen::Matrix<double, 1, Eigen::Dynamic> focalLengths_m = solutions.at(0);
	Eigen::Matrix<double, 1, Eigen::Dynamic> zb = solutions.at(1);
	Eigen::Matrix<double, 1, Eigen::Dynamic> zc = solutions.at(2);
	Eigen::Matrix<double, 1, Eigen::Dynamic> zd = solutions.at(3);
	
	// Check of the solutions
	if (DBG_M) std::cout << "This is focalLenghts:" << focalLengths_m << "\n";
	if (DBG_M) std::cout << "This is zb:" << zb << "\n";
	if (DBG_M) std::cout << "This is zc:" << zc << "\n";
	if (DBG_M) std::cout << "This is zd:" << zd << "\n";
	
	// number of focal lenghts possible found with check
	double lcnt = focalLengths_m.size();
	if (DBG_M) std::cout << "\nlcnt = " << lcnt <<".\n";
	
	//Vector for the translation vectors to be used in every loop
	std::vector<Eigen::Matrix<double, 3, 1>> translationVectors;
	
	// now we loop and go through every possible solution 
	for (int i=0; i<lcnt; i++)//i=1:lcnt
	{
		// print the loop cycle
		if (DBG_M) std::cout << "i=" << i <<"\n";
		
		// create p3d points in a camera coordinate system (using depths)
		Eigen::Matrix<double,3,4> p3dc = Eigen::Matrix<double,3,4>::Zero(3,4);
		Eigen::Matrix<double,3,1> temp = Eigen::Matrix<double,3,1>::Zero(3,1);
		
		// Check current points
		if (DBG_M) std::cout << "This is zb(i)=  "<< zb(i) <<"\n\n";
		if (DBG_M) std::cout << "This is zc(i)=  "<< zc(i) <<"\n\n";
		if (DBG_M) std::cout << "This is zd(i)=  "<< zd(i) <<"\n\n";
		if (DBG_M) std::cout << "This is focalLengths_m(i)=  "<< focalLengths_m(i) <<"\n\n";
		if (DBG_M) std::cout << "This is m2d\n"<< m2D <<"\n\n";
		
		// p3dc gets filled and verified at each step
		temp.col(0) << m2D.col(0),
		focalLengths_m(i);
		p3dc.col(0) = 1* temp;
		if (DBG_M) std::cout << "This is p3dc 0\n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(1),
		focalLengths_m(i);
		p3dc.col(1) = zb(i)* temp;
		if (DBG_M) std::cout << "This is p3dc 1 \n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(2),
		focalLengths_m(i);
		p3dc.col(2) = zc(i)* temp;
		if (DBG_M) std::cout << "This is p3dc 2 \n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(3),
		focalLengths_m(i);
		p3dc.col(3) = zd(i)* temp;
		if (DBG_M) std::cout << "This is p3dc 3 \n"<<p3dc <<"\n\n";
		
		// d matrix prepared and then calculated
		Eigen::Matrix<double,6,1> d_matrix = Eigen::Matrix<double,6,1>::Zero(6,1);
		
		// fix scale (recover 'za')
		double temp_double = 0;
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(1)) ).colwise().sum()(0);
		d_matrix(0) = sqrt(glab / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(2)) ).colwise().sum()(0);
		d_matrix(1) = sqrt(glac / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(2) = sqrt(glad / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(1) - p3dc.col(2)) ).colwise().sum()(0);
		d_matrix(3) = sqrt(glbc / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(1) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(4) = sqrt(glbd / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(2) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(5) = sqrt(glcd / temp_double);
		
		// d matrix checked
		if (DBG_M) std::cout << "This is d:\n"<< d_matrix<<"\n\n";
		
		//Check of d matrix
		double gta = d_matrix.colwise().sum()(0) / 6;
		if (DBG_M) std::cout << "This is gta= "<< gta<<"\n";
		p3dc = gta * p3dc;
		if (DBG_M) std::cout << "This is p3dc FINAL\n\n"<< p3dc <<"\n\n";
		
		//Vector for solution output		
		std::vector<Eigen::Matrix<double, 3,3>> rigidTrasform_results;
		
		// function is called and checked
		if (DBG_M) std::cout << "\n********************************************** BEGIN GETRIGIDTRANSFORM2**********************************************\n";
		int errorCheck = getRigidTransform2(M3D, p3dc, false, &rigidTrasform_results);
		if(errorCheck == -1)
		{
			return -1;
		}
		if (DBG_M) std::cout << "\n********************************************** END GETRIGIDTRANSFORM2 **********************************************\n";
		if (DBG_M) std::cout << "errorCheck is= "<< errorCheck << "\n\n";
		
		// Create temporary matrices to stuff them into the result vector
		Eigen::Matrix<double,3,3> Rr;
		Eigen::Matrix<double,3,3> tt_big;
		Eigen::Matrix<double,3,1> tt;
		Rr = rigidTrasform_results.at(0);
		tt_big = rigidTrasform_results.at(1);
		
		// Unpacking of translation vector 
		for (int i=0; i<3;i++)
		{
			tt(i) = tt_big(i,0);
		}
		
		if (DBG_M) std::cout << "RR is:\n"<< Rr << "\n\n";
		if (DBG_M) std::cout << "tt is:\n"<< tt << "\n\n";
		if (DBG_M) std::cout << "mean3d is:\n"<< mean3d << "\n\n";
		if (DBG_M) std::cout << "var is:\n"<< var << "\n\n";
		if (DBG_M) std::cout << "var*tt:\n"<< var*tt << "\n\n";
		if (DBG_M) std::cout << "Rr*mean3d:\n"<< Rr*mean3d << "\n\n";
		
		// Pushing out solutions
		rotationMatrices->push_back(Rr);
		temp = var*tt-Rr*mean3d;
		translationResult->push_back(temp);
		double focalToReturn = var2d*focalLengths_m(i);
		focalLengths_mResults->push_back(focalToReturn);
		
		// Print the solutions before sending them back
		if (DBG_M) std::cout << "focalToReturn:\n"<< focalToReturn << "\n";
		if (DBG_M) std::cout << "RR is:\n"<< Rr << "\n\n";
		if (DBG_M) std::cout << "tt is:\n"<< tt << "\n\n";
	}
	// Saying we are returning
	if (DBG_M) std::cout << "\n**************** END OF P4PF_M ****************\n";
	return 0;
}
/*********************************************/
/************ END second funciton ************/
/*********************************************/

/*********************************************/
/*********** BEGIN third funciton ************/
/*********************************************/
/* COMMENTS FROM THE MATLAB CODE
% P4P + unknown focal length
% given a set of 4x 2D<->3D correspondences, calculate camera pose and the
% camera focal length.
%
% by Martin Bujnak, (c)apr2008
%
%
% Please refer to the following paper, when using this code :
%
%      Bujnak, M., Kukelova, Z., and Pajdla, T. A general solution to the p4p
%      problem for camera with unknown focal length. CVPR 2008, Anchorage,
%      Alaska, USA, June 2008
%
%
% function [f R t] = P4Pf(m2D, M3D)
%
% input:
%
%  m2D - 2x4 matrix with 4x2D measuremets
%  M3D - 3x4 matrix with corresponding 3D points
%
% output:
%
%  f - vector with N focal lengths
%  R - 3x3xN matrix with N rotation matrices
%  t - 3xN matrix with N translation matrices
%
%  following equation holds for each solution
%
%      lambda*m2D = diag([f(i) f(i) 1])*[R(:,:,i) t(:,i)] * M3D
*/
/* 
INPUT: 	- m2D, 2x4 matrix
- M3D, 3x4 matrix
- a vector of double for returning focal lenghts
- a vector of matrices 3x3 for the rotation matrices
- a vector of Eigen vectors 3x1 for the translation vectors
OUTPUT:	- The given vectors filled
- the int return is used as error check: if 0 correct, if -1 an error was found
NOTE: This is similar to P4pf_M but with improved efficiency.
There are comments only on parts which differs from previous code.
*/ 
int P4pf::P4Pf(Eigen::Matrix<double, 2, 4> m2D,	Eigen::Matrix<double, 3, 4> M3D, std::vector<double>* focalLengthsResults, std::vector<Eigen::Matrix<double, 3,3>>* rotationMatrices, std::vector<Eigen::Matrix<double, 3,1>>* translationResult )
{
	/* SAME CODE OF P4PF_M*/
	if (DBG) std::cout << "\n**************** ENTERING P4PF ****************\n";
	if (DBG) std::cout << "this is M3D:\n" << M3D << "\n\n";
	if (DBG) std::cout << "this is M3D.transpose():\n" << M3D.transpose() << "\n\n";									//matrix trasposed but it's the same as before
	if (DBG) std::cout << "this is M3D.transpose().colwise().sum():\n" << M3D.transpose().colwise().sum() << "\n\n";
	
	Eigen::Matrix<double, 1, 3> mean3d = M3D.transpose().colwise().sum() /4;
	
	if (DBG) std::cout << "this is mean3d:\n" << mean3d << "\n\n";
	if (DBG) std::cout << "this is mean3d.replicate(1, 4):\n" << mean3d.replicate(1, 4) << "\n\n";
	if (DBG) std::cout << "\nmean3d created\n";
	if (DBG) std::cout << "this is mean3d:\n" << mean3d << "\n\n";
	if (DBG) std::cout << "this is mean3d.transpose():\n" << mean3d.transpose() << "\n\n";
	
	M3D = M3D - mean3d.transpose().replicate(1, 4);
	Eigen::Matrix<double, 3, 4> M3Dbis = Eigen::Matrix<double, 3, 4>::Zero(3,4);
	M3Dbis = this->pointWisePower(M3D);
	
	if (DBG) std::cout << "\nThis is M3D\n" << M3D << "\n";
	if (DBG) std::cout << "\nThis is M3D^2\n" << M3Dbis << "\n";
	if (DBG) std::cout << "\nThis is M3D colwise \n" << M3Dbis.colwise().sum() << "\n";
	if (DBG) std::cout << "\nThis is M3D colwise + sqrt \n" << M3Dbis.colwise().sum().cwiseSqrt() << "\n";
	if (DBG) std::cout << "\nThis is M3D colwise + sqrt + rowise sum\n" << M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum() << "\n";
	if (DBG) std::cout << "\nThis is M3D colwise + sqrt + rowise sum at 0\n" << M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum()(0) << "\n";
	
	double var = M3Dbis.colwise().sum().cwiseSqrt().rowwise().sum()(0);
	var = var/4;
	M3D = (1/var)*M3D;
	
	if (DBG) std::cout << "\nThis is var = " << var << "\n";
	if (DBG) std::cout << "\nThis is M3D final, after it was elaborated with var\n" << M3D << "\n";
	
	// scale 2D data
	double var2d = this->pointWisePower(m2D).colwise().sum().cwiseSqrt().rowwise().sum()(0);
	var2d = var2d/4;
	m2D = (1/var2d)*m2D;
	
	if (DBG) std::cout << "\n var2d extracted from m2D = " << var2d;
	if (DBG) std::cout << "\nThis is m2D final, after it was elaborated with var2D\n" << m2D << "\n";
	
	//coefficients of 5 reduced polynomials in these monomials mon
	double glab = this->pointWisePower( (M3D.col(0) - M3D.col(1)) ).colwise().sum()(0);
	double glac = this->pointWisePower( (M3D.col(0) - M3D.col(2)) ).colwise().sum()(0);
	double glad = this->pointWisePower( (M3D.col(0) - M3D.col(3)) ).colwise().sum()(0);
	double glbc = this->pointWisePower( (M3D.col(1) - M3D.col(2)) ).colwise().sum()(0);
	double glbd = this->pointWisePower( (M3D.col(1) - M3D.col(3)) ).colwise().sum()(0);
	double glcd = this->pointWisePower( (M3D.col(2) - M3D.col(3)) ).colwise().sum()(0);
	
	double	tol = 2.2204e-10;
	if (glbc*glbd*glcd*glab*glac*glad < tol)
	{
		return -1;
	}
	if (DBG) std::cout << "\n Passed the coefficient check, no zeros here\n";
	
	// Here the code icludes the C file from the paper.
	// Data needs to be prepared accordingly
	// Creation of a array of gl coefficients
	double gl_coef[6];
	gl_coef[0] = glab;
	gl_coef[1] = glac;
	gl_coef[2] = glad;
	gl_coef[3] = glbc;
	gl_coef[4] = glbd;
	gl_coef[5] = glcd;
	
	//other arrays of coefficients
	double a_coef[] = {(double) m2D(0,0),(double) m2D(1,0)};
	double b_coef[] = {(double) m2D(0,1),(double) m2D(1,1)};
	double c_coef[] = {(double) m2D(0,2), (double)m2D(1,2)};
	double d_coef[] = {(double) m2D(0,3),(double) m2D(1,3)};
	
	//Code of p4pfcode is used so we prepare matrices and solutions before.
	Eigen::Matrix<std::complex<double>, 10, 10> A = (Eigen::Matrix<std::complex<double>, 10, 10>::Zero(10,10));
	std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> solutions;
	if (DBG_M) std::cout << "solutions.size()= " << solutions.size() << "\n";
	
	//additional code is called and checked.
	this->P4pf::p4pfmex(gl_coef, a_coef, b_coef, c_coef, d_coef, &A);
	if (DBG) std::cout << "This is A\n\n" << A << "\n\n";	
	
	//matrix comes out transposed, so we traspose it back and use the code of p4pfcode
	Eigen::Matrix<std::complex<double>, 10, 10> temporaryA;
	temporaryA = A.transpose();
	A = temporaryA;
	// result check.
	if (DBG) std::cout << "This is A transposed\n\n" << A << "\n\n";
	
	/****************************** BEGIN CODE OF P4PCODE FUNCTION ************************/
	Eigen::Matrix<std::complex<double>, 1, 10> D_diag;
	Eigen::Matrix<std::complex<double>, 10, 10> D = (Eigen::Matrix<std::complex<double>, 10, 10>::Zero(10,10));
	Eigen::Matrix<std::complex<double>, 10, 10> V;
	Eigen::Matrix<std::complex<double>, 4, 10> sol;
	
	if (DBG) std::cout << "Ready for eigensolver\n";
	
	Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>, 10, 10>> es(A);
	
	D_diag =es.eigenvalues();
	if (DBG) std::cout << "This is D_diag\n\n" << D_diag << "\n\n";
	
	for (int i = 0; i < 10; i++)
	{
		D(i,i) = D_diag(i);
	}
	if (DBG) std::cout << "This is D\n\n" << D << "\n\n";
	
	V = es.eigenvectors();
	if (DBG) std::cout << "This is V\n" << V<< "\n\n";
	if (DBG) std::cout << "V.rows()" << V.rows() << "V.cols() " << V.cols() <<" V.size() " <<V.size() <<"\n";
	
	for(int i = 1; i<5; i++)
	{
		for(int j=0; j<10; j++)
		{
			sol(i-1,j) = V(i,j)/V(0,j);
			if (!(sol.real()(i-1,j) >= 0 || sol.real()(i-1,j) <0))
			{
				return -1;
			}
		}
	}
	std::vector<int> indexes;
	if (DBG) std::cout << "Now Entering the index search loop\n";
	
	for(int j = 0; j<10; j++)
	{
		if (DBG) std::cout << ". j = " << j <<"\n";
		if (DBG) std::cout << "real = " << sol.real()(3,j) <<"\n";
		if (DBG) std::cout << "imag = " << sol.imag()(3,j) <<"\n";
		if (DBG) std::cout << "imaginary check = " << (sol.imag()(3,j)<0.00001 && sol.imag()(3,j) > -0.00001) <<"\n";
		
		if ((sol.imag()(3,j)<0.00001 && sol.imag()(3,j) > -0.00001) ) //solution is non imaginary
		{
			if (DBG) std::cout << "positive check = " << (sol.real()(3,j)>0.00001) <<"\n";
			if (sol.real()(3,j)>0.00001)
			{
				if (DBG) std::cout << "I got this\n";
				indexes.push_back(j);
			}
		}
	}
	
	if (DBG) std::cout << "This is sol\n" << sol<< "\n";
	if (DBG) std::cout << "Found focal lengths, now saving solutions\n";
	
	Eigen::Matrix<double, 1, Eigen::Dynamic> f;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zd;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zc;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zb;
	
	f.resize(1,indexes.size());
	zd.resize(1,indexes.size());
	zc.resize(1,indexes.size());
	zb.resize(1,indexes.size());
	
	if (DBG) std::cout << "sol.rows() = " << sol.rows() << ". sol.cols() = " << sol.cols() <<". sol.size() = " <<sol.size() <<".\n";
	
	if (DBG) std::cout << "indexes.size() = " << indexes.size() << "\n\n";
	
	//go to positions marked and take results
	for (int i = 0; i < indexes.size(); i++)
	{
		if (DBG) std::cout << "i = " << i<< "\n\n";
		
		f(i)  = sqrt(sol.real()(3, indexes.at(i)));
		zd(i) = sol.real()(0, indexes.at(i));
		zc(i) = sol.real()(1, indexes.at(i));
		zb(i) = sol.real()(2, indexes.at(i));
		if (DBG) std::cout << "f(" <<i<< ")= " << f(i) <<"\n";
		if (DBG) std::cout << "zd(" <<i<< ")= " << zd(i) <<"\n";
		if (DBG) std::cout << "zc(" <<i<< ")= " << zc(i) <<"\n";
		if (DBG) std::cout << "zb(" <<i<< ")= " << zb(i) <<"\n";
	}
	
	if (DBG) std::cout << "Done. Now Onto the next part\n";
	/****************************** END CODE OF P4PCODE FUNCTION ************************/
	
	// compatibility with previous code, we change name
	Eigen::Matrix<double, 1, Eigen::Dynamic> focalLengths = f;
	if (DBG) std::cout << "This is focalLenghts:\n" << focalLengths <<"\n";
	
	double lcnt = focalLengths.size();
	std::vector<Eigen::Matrix<double, 3, 1>> translationVectors;
	if (DBG) std::cout << "\nlcnt = " << lcnt <<".\n";
	
	for (int i=0; i<lcnt; i++)
	{
		if (DBG) std::cout << "i=" << i <<"\n";
		
		Eigen::Matrix<double,3,4> p3dc = Eigen::Matrix<double,3,4>::Zero(3,4);
		Eigen::Matrix<double,3,1> temp = Eigen::Matrix<double,3,1>::Zero(3,1);
		
		if (DBG) std::cout << "This is zb(i)=  "<< zb(i) <<"\n\n";
		if (DBG) std::cout << "This is zc(i)=  "<< zc(i) <<"\n\n";
		if (DBG) std::cout << "This is zd(i)=  "<< zd(i) <<"\n\n";
		if (DBG) std::cout << "This is focalLengths(i)=  "<< focalLengths(i) <<"\n\n";
		if (DBG) std::cout << "This is m2d\n"<< m2D <<"\n\n";
		
		temp.col(0) << m2D.col(0),
		focalLengths(i);
		p3dc.col(0) = 1* temp;
		if (DBG) std::cout << "This is p3dc 0\n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(1),
		focalLengths(i);
		p3dc.col(1) = zb(i)* temp;
		if (DBG) std::cout << "This is p3dc 1 \n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(2),
		focalLengths(i);
		p3dc.col(2) = zc(i)* temp;
		if (DBG) std::cout << "This is p3dc 2 \n"<<p3dc <<"\n\n";
		
		temp.col(0) << m2D.col(3),
		focalLengths(i);
		p3dc.col(3) = zd(i)* temp;
		
		if (DBG) std::cout << "This is p3dc 3 \n"<<p3dc <<"\n\n";
		if (DBG) std::cout << "temp and p3dc full\n\n";
		Eigen::Matrix<double,6,1> d_matrix = Eigen::Matrix<double,6,1>::Zero(6,1);
		
		if (DBG) std::cout << "This is temp\n"<< temp <<"\n\n";
		
		double temp_double = 0;
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(1)) ).colwise().sum()(0);
		d_matrix(0) = sqrt(glab / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(2)) ).colwise().sum()(0);
		d_matrix(1) = sqrt(glac / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(0) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(2) = sqrt(glad / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(1) - p3dc.col(2)) ).colwise().sum()(0);
		d_matrix(3) = sqrt(glbc / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(1) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(4) = sqrt(glbd / temp_double);
		
		temp_double = this->pointWisePower( (p3dc.col(2) - p3dc.col(3)) ).colwise().sum()(0);
		d_matrix(5) = sqrt(glcd / temp_double);
		
		double gta = d_matrix.colwise().sum()(0) / 6;
		p3dc = gta * p3dc;
		
		if (DBG) std::cout << "This is d:\n"<< d_matrix<<"\n\n";
		if (DBG) std::cout << "This is gta= "<< gta<<"\n";
		if (DBG) std::cout << "This is p3dc FINAL\n\n"<< p3dc <<"\n\n";
		
		std::vector<Eigen::Matrix<double, 3,3>> rigidTrasform_results;
		
		if (DBG) std::cout << "\n********************************************** BEGIN GETRIGIDTRANSFORM2**********************************************\n";
		
		int errorCheck = getRigidTransform2(M3D, p3dc, false, &rigidTrasform_results);
		if(errorCheck == -1)
		{
			return -1;
		}
		if (DBG) std::cout << "\n********************************************** END GETRIGIDTRANSFORM2 **********************************************\n";
		if (DBG) std::cout << "errorCheck is= "<< errorCheck << "\n\n";
		
		Eigen::Matrix<double,3,3> Rr;
		Eigen::Matrix<double,3,3> tt_big;
		Eigen::Matrix<double,3,1> tt;
		Rr = rigidTrasform_results.at(0);
		tt_big = rigidTrasform_results.at(1);
		
		for (int i=0; i<3;i++)
		{
			tt(i) = tt_big(i,0);
		}
		
		if (DBG) std::cout << "RR is:\n"<< Rr << "\n\n";
		if (DBG) std::cout << "tt is:\n"<< tt << "\n\n";
		if (DBG) std::cout << "mean3d is:\n"<< mean3d << "\n\n";
		if (DBG) std::cout << "var is:\n"<< var << "\n\n";
		if (DBG) std::cout << "var*tt:\n"<< var*tt << "\n\n";
		if (DBG) std::cout << "Rr*mean3d:\n"<< Rr*(mean3d.transpose()) << "\n\n";
		
		rotationMatrices->push_back(Rr);
		temp = var*tt-Rr*(mean3d.transpose());
		translationResult->push_back(temp);
		double focalToReturn = var2d*focalLengths(i);
		focalLengthsResults->push_back(focalToReturn);
		
		if (DBG) std::cout << "focalToReturn:\n"<< focalToReturn << "\n";
		if (DBG) std::cout << "RR is:\n"<< Rr << "\n\n";
		if (DBG) std::cout << "tt is:\n"<< tt << "\n\n";
		
	}
	if (DBG) std::cout << "\n**************** END OF P4PF_M ****************\n";
	return 0;
}
/*********************************************/
/************* END third funciton ************/
/*********************************************/

/* FUNCTION THAT FILLS A GIVEN MATRIX IN PARTICULAR SPARSE COEFFICIENTS 

INPUT: - the 78x88 matrix to fill, as a pointer
- the array of coefficients
- the value to store inside the matrix cells
- the size of the array passed
OUTPUT: 0 when done
*/

int P4pf::fillMatrix(Eigen::Matrix<double,78,88> *matrix, int *arr, double value, int size)
{
	// check which size was received
	if (DBG_FILL) std::cout << "size received = " << size << "\n";
	for ( int i = 0; i < size; i++)
	{
		// analyse every cycle
		if (DBG_FILL) std::cout << "i = " << i << "\n";
		if (DBG_FILL) std::cout << "arr[i] = " << arr[i] << "\n";
		if (DBG_FILL) std::cout << "value = " << value << "\n";
		
		//column and rows calculated since indexes inside the arrays are in matlab format
		// so they start from 1.
		int col = (arr[i]-1) / 78;
		int row = (arr[i] -1 ) % 78;
		
		(*matrix)(row,col)= value;
		
		// Double check if it was inserted correctly, you never know
		if (DBG_FILL) std::cout << "Matrix value inserted = " << (*matrix)(row,col) << "\n\n";
	}
	return 0;
}

/*********************************************/
/*********** BEGIN fourth funciton ************/
/*********************************************/

/* MATLAB COMMENTS
% P4P + unknown focal length helper
% given a set of 4x 2D<->3D correspondences, calculate focal length and 3 unknown depths
%
% by Martin Bujnak, (c)june2011
%
%
% Please refer to the following paper, when using this code :
%
%      Bujnak, M., Kukelova, Z., and Pajdla, T. A general solution to the p4p
%      problem for camera with unknown focal length. CVPR 2008, Anchorage,
%      Alaska, USA, June 2008
%
% [f zb zc zd] = p4pfcode(glab, glac, glad, glbc, glbd, glcd, a1, a2, b1, b2, c1, c2, d1, d2)
%
% input:
%   glXY - ||X-Y||^2 - quadratic distances between 3D points X and Y
%   a1 (a2) = x (resp y) measurement of the first 2D point
%   b1 (b2) = x (resp y) measurement of the second 2D point
%   c1 (c2) = x (resp y) measurement of the third 2D point
%   d1 (d2) = x (resp y) measurement of the fourth 2D point
%
% output:
%   f - estimated focal lengths
%   zb, zc, zd  - depths of points b, c, d (depth of 'a' is fixed to 1)
%
*/
/* FUNCTION THAT FILLS A GIVEN MATRIX IN PARTICULAR SPARSE COEFFICIENTS 

INPUT: - values needed for the calculations
- vector for the solutions
OUTPUT: - 0 if executed correctly
- -1 if errors
*/
int P4pf::p4pfcode(double glab, double glac, double glad,
	double glbc, double glbd, double glcd,
	double a1, double a2, double b1, double b2,
	double c1, double c2, double d1, double d2, std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>>* solutions)
{
	//Check if variables were received correctly
	if (DBG_P4PFCODE) std::cout << "\n**************** ENTERED P4PFCODE ****************\n";
	if (DBG_P4PFCODE) std::cout << "\n**************** PARAMETERS:\n";
	if (DBG_P4PFCODE) std::cout << "glab = " << glab << "\n";
	if (DBG_P4PFCODE) std::cout << "glac = " << glac << "\n";
	if (DBG_P4PFCODE) std::cout << "glad = " << glad << "\n";
	if (DBG_P4PFCODE) std::cout << "glbc = " << glbc << "\n";
	if (DBG_P4PFCODE) std::cout << "glbd = " << glbd << "\n";
	if (DBG_P4PFCODE) std::cout << "glcd = " << glcd << "\n";
	if (DBG_P4PFCODE) std::cout << "a1 = " << a1 << "\n";
	if (DBG_P4PFCODE) std::cout << "a2 = " << a2 << "\n";
	if (DBG_P4PFCODE) std::cout << "b1 = " << b1 << "\n";
	if (DBG_P4PFCODE) std::cout << "b2 = " << b2 << "\n";
	if (DBG_P4PFCODE) std::cout << "c1 = " << c1 << "\n";
	if (DBG_P4PFCODE) std::cout << "c2 = " << c2 << "\n";
	if (DBG_P4PFCODE) std::cout << "d1 = " << d1 << "\n";
	if (DBG_P4PFCODE) std::cout << "d2 = " << d2 << "\n";
	if (DBG_P4PFCODE) std::cout << "\n**************** NOW THE CALCULATIONS ****************\n";
	
	// precalculate polynomial equations coefficients
	Eigen::Matrix<double, 78,88> M = Eigen::Matrix<double,78,88>::Zero(78,88);
	if (DBG_P4PFCODE) std::cout << "Start filling matrix\n";
	
	//The matrix gets filled
	int id_01[] = {72, 149, 520, 597, 752, 829, 1062, 1217, 1528, 1895, 2050, 2127, 2360, 2515, 2904, 3439, 3594, 3983, 4993};  //size 19
	fillMatrix(&M, id_01, 1, sizeof(id_01)/sizeof(id_01[0]));
	int id_02[] = {384, 461, 988, 1299, 1454, 1609, 1686, 1841, 2830, 2985, 3140, 3295, 3372, 4219, 4374, 4607, 5539}; //17
	fillMatrix(&M, id_02,0.5/glad*glbc-0.5*glab/glad-0.5*glac/glad, sizeof(id_02 )/sizeof(id_01[0]));
	int id_03[] = {618, 929, 1924, 2235, 2390, 2545, 2778, 2933, 3244, 3454, 3609, 3842, 3997, 4308, 4843, 4998, 5231, 5851  }; //18
	fillMatrix(&M, id_03,-1, sizeof(id_03 )/sizeof(id_01[0]));
	int id_04[] = {696, 1007, 2002, 2313, 2468, 2623, 2856, 3011, 3322, 3532, 3687, 3920, 4075, 4386, 4921, 5076, 5309, 5929  };//18
	fillMatrix(&M, id_04, -1, sizeof(id_04 )/sizeof(id_01[0]));
	int id_05[] = {774, 1085, 2080, 2391, 2546, 2701, 2934, 3089, 3455, 3610, 3765, 3998, 4153, 4464, 4999, 5154, 5387, 6007 }; //18
	fillMatrix(&M, id_05,c2*b2+c1*b1, sizeof(id_05 )/sizeof(id_01[0]));
	int id_06[] = {1008, 1319, 2314, 2858, 3013, 3168, 3323, 3400, 3922, 4077, 4232, 4387, 4620, 5310, 5543, 6163  }; //18
	fillMatrix(&M, id_06, glac/glad-1/glad*glbc+glab/glad, sizeof(id_06 )/sizeof(id_01[0]));
	int id_07[] = {1476, 1709, 2860, 3171, 3326, 3402, 4235, 4390, 4545, 4622, 4777, 5545, 5700, 5777, 6397  }; //15
	fillMatrix(&M, id_07, 0.5/glad*glbc*pow(d2,2)-0.5*glab/glad*pow(d2,2)-0.5*glac/glad*pow(d2,2)-0.5*glac/glad*pow(d1,2)+0.5/glad*glbc*pow(d1,2)-0.5*glab/glad*pow(d1,2), sizeof(id_07 )/sizeof(id_01[0]));
	int id_08[] = {2334, 3950, 4105, 4260, 4415, 4648, 4936, 5091, 5323, 5556, 5934, 6167, 6475}; //13
	fillMatrix(&M, id_08, 1-0.5*glac/glad-0.5*glab/glad+0.5/glad*glbc, sizeof(id_08 )/sizeof(id_01[0]));
	int id_09[] = {2412, 2801, 3484, 3873, 4028, 4183, 4338, 4493, 4726, 4859, 5014, 5169, 5246, 5401, 5634, 5857, 6012, 6245, 6553 }; //19
	fillMatrix(&M, id_09, -b1*a1-a2*b2, sizeof(id_09 )/sizeof(id_01[0]));
	int id_10[] = {2490, 2879, 3562, 3951, 4106, 4416, 4571, 4804, 4937, 5092, 5324, 5479, 5712, 5935, 6090, 6323, 6631  }; //17
	fillMatrix(&M, id_10, -c2*a2-c1*a1, sizeof(id_10 )/sizeof(id_01[0]));
	int id_11[] = {2880, 3191, 3952, 4263, 4418, 4573, 4650, 4805, 5326, 5481, 5558, 5713, 5790, 6169, 6324, 6401, 6709  }; //17
	fillMatrix(&M, id_11, -a1/glad*glbc*d1+a1*glac/glad*d1+glac/glad*a2*d2+a1*glab/glad*d1-1/glad*glbc*a2*d2+glab/glad*a2*d2, sizeof(id_11 )/sizeof(id_01[0]));
	int id_12[] = {3972, 4283, 4966, 5354, 5509, 5586, 5741, 5818, 5950, 6105, 6182, 6337, 6414, 6481, 6636, 6713, 6787 }; //17
	fillMatrix(&M, id_12, pow(a2,2)+pow(a1,2)-0.5*glac/glad*pow(a2,2)-0.5*pow(a1,2)*glac/glad+0.5/glad*glbc*pow(a2,2)-0.5*pow(a1,2)*glab/glad+0.5*pow(a1,2)/glad*glbc-0.5*glab/glad*pow(a2,2), sizeof(id_12 )/sizeof(id_01[0]));
	int id_13[] = {74, 229, 527, 682, 759, 836, 1147, 1224, 1613, 1980, 2057, 2134, 2445, 2522, 2599, 2988, 3520, 3597, 4064, 5072  }; //20
	fillMatrix(&M, id_13, 1, sizeof(id_13 )/sizeof(id_01[0]));
	int id_14[] = {308, 463, 917, 1306, 1383, 1538, 1693, 1770, 2759, 2914, 3147, 3224, 3301, 3378, 4222, 4299, 4610, 5540  }; //18
	fillMatrix(&M, id_14, -glac/glad, sizeof(id_14 )/sizeof(id_01[0]));
	int id_15[] = {620, 1009, 1931, 2320, 2397, 2552, 2863, 2940, 3329, 3461, 3616, 3927, 4004, 4081, 4392, 4924, 5001, 5312, 5930  };//19
	fillMatrix(&M, id_15, -2, sizeof(id_15 )/sizeof(id_01[0]));
	int id_16[] = {776, 1165, 2087, 2476, 2553, 2708, 3019, 3096, 3540, 3617, 3772, 4083, 4160, 4548, 5080, 5157, 5468, 6086  }; //18
	fillMatrix(&M, id_16, pow(c1,2)+pow(c2,2), sizeof(id_16 )/sizeof(id_01[0]));
	int id_17[] = {932, 1321, 2243, 2787, 2942, 3175, 3252, 3407, 3851, 4006, 4239, 4316, 4393, 4626, 5235, 5546, 6164  }; //17
	fillMatrix(&M, id_17, 2*glac/glad, sizeof(id_17 )/sizeof(id_01[0]));
	int id_18[] = {1400, 1711, 2789, 3178, 3255, 3409, 4242, 4319, 4474, 4629, 4706, 4783, 5548, 5625, 5780, 6398  };
	fillMatrix(&M, id_18, -glac/glad*pow(d1,2)-glac/glad*pow(d2,2), sizeof(id_18 )/sizeof(id_01[0]));
	int id_19[] = {2258, 3879, 4034, 4267, 4344, 4655, 4865, 5020, 5252, 5329, 5562, 5859, 6170, 6476  };
	fillMatrix(&M, id_19,-glac/glad+1, sizeof(id_19 )/sizeof(id_01[0]));
	int id_20[] = {2414, 2881, 3491, 3958, 4035, 4190, 4423, 4500, 4811, 4944, 5021, 5176, 5331, 5408, 5485, 5718, 5938, 6015, 6326, 6632  };
	fillMatrix(&M, id_20,-2*c2*a2-2*c1*a1, sizeof(id_20 )/sizeof(id_01[0]));
	int id_21[] = {2804, 3193, 3881, 4270, 4347, 4502, 4657, 4734, 5255, 5410, 5565, 5642, 5719, 5796, 6172, 6249, 6404, 6710  };
	fillMatrix(&M, id_21,2*a1*glac/glad*d1+2*glac/glad*a2*d2, sizeof(id_21)/sizeof(id_01[0]));
	int id_22[] = {3896, 4285, 4895, 5283, 5438, 5593, 5670, 5825, 5879, 6034, 6189, 6266, 6343, 6420, 6484, 6561, 6716, 6788  };
	fillMatrix(&M, id_22,-glac/glad*pow(a2,2)+pow(a2,2)+pow(a1,2)-pow(a1,2)*glac/glad, sizeof(id_22 )/sizeof(id_01[0]));
	int id_23[] = {154, 309, 609, 920, 1075, 1386, 2220, 2375, 2530, 2763, 2918, 3229, 3835, 3990, 4301, 5229  };
	fillMatrix(&M, id_23,1, sizeof(id_23 )/sizeof(id_01[0]));
	int id_24[] = {388, 465, 999, 1310, 1465, 1698, 2843, 2998, 3153, 3308, 3385, 4225, 4380, 4613, 5541  };
	fillMatrix(&M, id_24,0.5/glad*glbd-0.5-0.5*glab/glad, sizeof(id_24 )/sizeof(id_01[0]));
	int id_25[] = {622, 933, 1935, 2246, 2401, 2790, 3467, 3622, 3855, 4010, 4321, 4849, 5004, 5237, 5853  };
	fillMatrix(&M, id_25, -1, sizeof(id_25 )/sizeof(id_01[0]));
	int id_26[] = {1012, 1323, 2325, 2869, 3180, 3935, 4090, 4245, 4400, 4633, 5316, 5549, 6165  };
	fillMatrix(&M, id_26, glab/glad-1/glad*glbd, sizeof(id_26 )/sizeof(id_01[0]));
	int id_27[] = {1090, 1401, 2403, 2792, 2947, 3258, 3858, 4013, 4168, 4323, 4478, 4711, 5239, 5394, 5627, 6243  };
	fillMatrix(&M, id_27,d2*b2+b1*d1, sizeof(id_27 )/sizeof(id_01[0]));
	int id_28[] = {1480, 1713, 2871, 3182, 3337, 3414, 4248, 4403, 4558, 4635, 4790, 5551, 5706, 5783, 6399  };
	fillMatrix(&M, id_28,-0.5*glab/glad*pow(d2,2)-0.5*glab/glad*pow(d1,2)+0.5/glad*glbd*pow(d2,2)+0.5/glad*glbd*pow(d1,2)-0.5*pow(d2,2)-0.5*pow(d1,2), sizeof(id_28 )/sizeof(id_01[0]));
	int id_29[] = {2338, 3961, 4272, 4949, 5104, 5336, 5569, 5940, 6173, 6477  };
	fillMatrix(&M, id_29,-0.5*glab/glad+0.5/glad*glbd+0.5, sizeof(id_29 )/sizeof(id_01[0]));
	int id_30[] = {2416, 2805, 3495, 3884, 4039, 4350, 4872, 5027, 5182, 5259, 5414, 5647, 5863, 6018, 6251, 6555  };
	fillMatrix(&M, id_30, -a2*b2-b1*a1, sizeof(id_30 )/sizeof(id_01[0]));
	int id_31[] = {2884, 3195, 3963, 4274, 4429, 4662, 5339, 5494, 5571, 5726, 5803, 6175, 6330, 6407, 6711  };
	fillMatrix(&M, id_31, -a1/glad*glbd*d1+a1*glab/glad*d1+glab/glad*a2*d2-1/glad*glbd*a2*d2, sizeof(id_31 )/sizeof(id_01[0]));
	int id_32[] = {3976, 4287, 4977, 5365, 5598, 5963, 6118, 6195, 6350, 6427, 6487, 6642, 6719, 6789  };
	fillMatrix(&M, id_32, 0.5/glad*glbd*pow(a2,2)+0.5*pow(a1,2)/glad*glbd-0.5*glab/glad*pow(a2,2)-0.5*pow(a1,2)*glab/glad+0.5*pow(a1,2)+0.5*pow(a2,2), sizeof(id_32 )/sizeof(id_01[0]));
	int id_33[] = {234, 389, 694, 849, 1004, 1159, 1470, 1547, 1624, 2307, 2384, 2461, 2538, 2615, 2848, 2925, 3002, 3313, 3917, 3994, 4071, 4382, 5308  };
	fillMatrix(&M, id_33, 1, sizeof(id_33 )/sizeof(id_01[0]));
	int id_34[] = {390, 467, 1006, 1239, 1316, 1471, 1704, 1781, 1858, 2774, 2851, 2928, 3005, 3160, 3237, 3314, 3391, 4229, 4306, 4383, 4616, 5542  };
	fillMatrix(&M, id_34,-0.5*glac/glad+0.5*glcd/glad-0.5, sizeof(id_34 )/sizeof(id_01[0]));
	int id_35[] = {702, 1013, 2020, 2175, 2330, 2485, 2874, 2951, 3028, 3476, 3553, 3630, 3707, 3940, 4017, 4094, 4405, 4931, 5008, 5085, 5318, 5932  };
	fillMatrix(&M, id_35,-1, sizeof(id_35 )/sizeof(id_01[0]));
	int id_36[] = {1014, 1325, 2332, 2565, 2875, 3186, 3263, 3340, 3866, 3943, 4020, 4097, 4252, 4329, 4406, 4639, 5242, 5319, 5552, 6166  };
	fillMatrix(&M, id_36,glac/glad-glcd/glad, sizeof(id_36 )/sizeof(id_01[0]));
	int id_37[] = {1170, 1481, 2488, 2721, 2876, 3031, 3342, 3945, 4022, 4099, 4176, 4408, 4485, 4562, 4795, 5321, 5398, 5475, 5708, 6322  };
	fillMatrix(&M, id_37,c1*d1+c2*d2, sizeof(id_37 )/sizeof(id_01[0]));
	int id_38[] = {1482, 1715, 2878, 3111, 3188, 3343, 3420, 4257, 4334, 4411, 4488, 4565, 4642, 4719, 4796, 5555, 5632, 5709, 5786, 6400  };
	fillMatrix(&M, id_38,0.5*glcd/glad*pow(d1,2)-0.5*glac/glad*pow(d2,2)-0.5*glac/glad*pow(d1,2)-0.5*pow(d1,2)-0.5*pow(d2,2)+0.5*glcd/glad*pow(d2,2), sizeof(id_38 )/sizeof(id_01[0]));
	int id_39[] = {2340, 3657, 3967, 4278, 4355, 4432, 4880, 4957, 5034, 5111, 5265, 5342, 5575, 5866, 5943, 6176, 6478  };
	fillMatrix(&M, id_39,-0.5*glac/glad+0.5+0.5*glcd/glad, sizeof(id_39 )/sizeof(id_01[0]));
	int id_40[] = {2496, 2885, 3580, 3813, 3968, 4123, 4434, 4511, 4588, 4959, 5036, 5113, 5190, 5344, 5421, 5498, 5731, 5945, 6022, 6099, 6332, 6634  };
	fillMatrix(&M, id_40,-c1*a1-c2*a2, sizeof(id_40 )/sizeof(id_01[0]));
	int id_41[] = {2886, 3197, 3970, 4203, 4280, 4435, 4668, 4745, 4822, 5270, 5347, 5424, 5501, 5578, 5655, 5732, 5809, 6179, 6256, 6333, 6410, 6712  };
	fillMatrix(&M, id_41,glac/glad*a2*d2+a1*glac/glad*d1-glcd/glad*a2*d2-glcd*a1/glad*d1, sizeof(id_41 )/sizeof(id_01[0]));
	int id_42[] = {3978, 4289, 4984, 5217, 5371, 5604, 5681, 5758, 5894, 5971, 6048, 6125, 6202, 6279, 6356, 6433, 6491, 6568, 6645, 6722, 6790  };
	fillMatrix(&M, id_42, 0.5*pow(a1,2)+0.5*pow(a2,2)-0.5*glac/glad*pow(a2,2)+0.5*glcd*pow(a1,2)/glad-0.5*pow(a1,2)*glac/glad+0.5*glcd/glad*pow(a2,2), sizeof(id_42 )/sizeof(id_01[0]));
	
	//Code for printing part of the matrix
	/*
		for (int i=0; i<78; i++)
		{
			std::cout << "row " << i << "\n";
			for (int j = 0; j<20; j++)
			{
				std::cout << M(i,j) << " ";
			}
			std::cout << "\n";
		}
	*/
	
	// Check M
	if (DBG_P4PFCODE) std::cout << "Matrix Filled\n";
	if (DBG_P4PFCODE) std::cout << "\n\nTHIS IS M\n" << M;
	if (DBG_P4PFCODE) std::cout << "\nCreating left and right matrices\n";
	
	//Creation of the two halfs
	Eigen::Matrix<double, 78, 78> M_left;
	Eigen::Matrix<double, 78, 10> M_right;
	if (DBG_P4PFCODE) std::cout << "Left and right created\n";
	for (int i = 0; i<78; i++)
	{
		for(int j = 0; j<78;j++)
		{
			M_left(i,j) = M(i,j);
		}
	}
	for (int i = 0; i<78; i++)
	{
		for(int j = 0; j<10;j++)
		{
			M_right(i,j) = M(i,j+78);
		}
	}
	
	if (DBG_P4PFCODE) std::cout << "Left and right matrices filled\n";
	
	Eigen::Matrix<double, 78, 10> Mr;
	Mr =M_left.lu().solve(M_right);  // Stable and fast. #include <Eigen/LU>
	//std::cout << "questa e M\n\n" << M << "\n\n";
	if (DBG_P4PFCODE) std::cout << "This is M_left\n\n" << M_left << "\n\n";
	if (DBG_P4PFCODE) std::cout << "This is M_right\n\n" << M_right << "\n\n";


	if (DBG_P4PFCODE) std::cout << "This is M_left\n\n" << M_left << "\n\n";
	if (DBG_P4PFCODE) std::cout << "This is M_right\n\n" << M_right << "\n\n";
	
	// Action matrix created and filled
	Eigen::Matrix<std::complex<double>, 10, 10> A = (Eigen::Matrix<std::complex<double>, 10, 10>::Zero(10,10));
	Eigen::Matrix<double, 1, 10> amcols;
	amcols << 9, 8, 7, 6, 5, 4, 3, 2, 1, 0;
	
	A(0, 1) = 1;
	A(1, 5) = 1;
	A(2, 6) = 1;
	A(3, 7) = 1;
	A(4, 8) = 1;	
	for(int i = 5; i < 10; i++)
	{
		for(int j = 0; j< 10;j++)
		{
			A(i,j) = -Mr(79-i, amcols(j));
		}
	}
	// Check A
	if (DBG_P4PFCODE) std::cout << "This is A\n\n" << A << "\n\n";
	
	//Now we need to solve and find eigenvector and eigenvalues
	// Used ComplexEigenSolver
	Eigen::Matrix<std::complex<double>, 1, 10> D_diag;
	Eigen::Matrix<std::complex<double>, 10, 10> D = (Eigen::Matrix<std::complex<double>, 10, 10>::Zero(10,10));
	Eigen::Matrix<std::complex<double>, 10, 10> V;
	Eigen::Matrix<std::complex<double>, 4, 10> sol;
	
	if (DBG_P4PFCODE) std::cout << "Ready for eigensolver\n";
	
	// Eigensolver used here
	Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>, 10, 10>> es(A);
	
	// Eigenvalues are only given back as a vector but we want them in a diagonal matrix
	D_diag =es.eigenvalues();
	if (DBG_P4PFCODE) std::cout << "This is D_diag\n\n" << D_diag << "\n\n";

	for (int i = 0; i < 10; i++)
	{
		D(i,i) = D_diag(i);
	}
	// Check result
	if (DBG_P4PFCODE) std::cout << "This is D_diag\n\n" << D_diag << "\n\n";
	if (DBG_P4PFCODE) std::cout << "This is D\n\n" << D << "\n\n";
	
	// We extract matrix of eigenvectors
	V = es.eigenvectors();
	if (DBG_P4PFCODE) std::cout << "This is V\n" << V<< "\n\n";
	if (DBG_P4PFCODE) std::cout << "V.rows()" << V.rows() << "V.cols() " << V.cols() <<" V.size() " <<V.size() <<"\n";
	
	//Personal method to check if it's Nan, isnan wasn't working with complex numbers. 
	for(int i = 1; i<5; i++)
	{
		for(int j=0; j<10; j++)
		{
			sol(i-1,j) = V(i,j)/V(0,j);
			if (!(sol.real()(i-1,j) >= 0 || sol.real()(i-1,j) <0))// personal method to check if NAN.
			{
				return -1;
			}	
		}
	}
	//OK, solution matrix does not hava NAN here.
	
	// We need to find REAL and compatible solution. 
	// To do this, we find indexes and save them because they are used later 
	std::vector<int> indexes;
	if (DBG_P4PFCODE) std::cout << "Now Entering the index search loop\n";
	for(int j = 0; j<10; j++)
	{
		//analyse 4th row of SOL and
		if (DBG_P4PFCODE) std::cout << ". j = " << j <<"\n";
		if (DBG_P4PFCODE) std::cout << "real = " << sol.real()(3,j) <<"\n";
		if (DBG_P4PFCODE) std::cout << "imag = " << sol.imag()(3,j) <<"\n";
		if (DBG_P4PFCODE) std::cout << "imaginary check = " << (sol.imag()(3,j)<0.00001 && sol.imag()(3,j) > -0.00001) <<"\n";
		
		// We check if a solution is immaginary we compare if is below a threshold.
		// This was tricky because the EigenSolver used to give 1e-6 values giving a non-diagonal matrix 
		if ((sol.imag()(3,j)<0.00001 && sol.imag()(3,j) > -0.00001) ) //solution is non imaginary
		{
			if (DBG_P4PFCODE) std::cout << "positive check = " << (sol.real()(3,j)>0.00001) <<"\n";
			if (sol.real()(3,j)>0.00001)	//Check if a solution is positive
			{
				//here soultions are positive and real, possible good solution.
				// We save it and print that we got it
				indexes.push_back(j);
				if (DBG_P4PFCODE) std::cout << "Got it\n";
			}
		}
	}
	
	//Found solution, we print it
	if (DBG_P4PFCODE) std::cout << "This is sol\n" << sol<< "\n";
	
	//here I know how many correct focal lengths I have
	//create matrices w/ correct length
	Eigen::Matrix<double, 1, Eigen::Dynamic> f;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zd;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zc;
	Eigen::Matrix<double, 1, Eigen::Dynamic> zb;
	
	// Dynamic matrices need to be sized correctly
	f.resize(1,indexes.size());
	zd.resize(1,indexes.size());
	zc.resize(1,indexes.size());
	zb.resize(1,indexes.size());
	
	// print which sizes are we talking about.
	if (DBG_P4PFCODE) std::cout << "sol.rows() = " << sol.rows() << ". sol.cols() = " << sol.cols() <<". sol.size() = " <<sol.size() <<".\n";
	if (DBG_P4PFCODE) std::cout << "indexes.size() = " << indexes.size() << "\n\n";
	
	//go to positions marked and take results
	for (int i = 0; i < indexes.size(); i++)
	{		
		f(i)  = sqrt(sol.real()(3, indexes.at(i)));
		zd(i) = sol.real()(0, indexes.at(i));
		zc(i) = sol.real()(1, indexes.at(i));
		zb(i) = sol.real()(2, indexes.at(i));
		
		if (DBG_P4PFCODE) std::cout << "i = " << i<< "\n\n";
		if (DBG_P4PFCODE) std::cout << "f(" <<i<< ")= " << f(i) <<"\n";
		if (DBG_P4PFCODE) std::cout << "zd(" <<i<< ")= " << zd(i) <<"\n";
		if (DBG_P4PFCODE) std::cout << "zc(" <<i<< ")= " << zc(i) <<"\n";
		if (DBG_P4PFCODE) std::cout << "zb(" <<i<< ")= " << zb(i) <<"\n";
	}
	
	//Passing solution back
	if (DBG_P4PFCODE) std::cout << "Done. Now passing solutions\n";
	solutions->push_back(f);
	solutions->push_back(zb);
	solutions->push_back(zc);
	solutions->push_back(zd);
	if (DBG_P4PFCODE) std::cout << "\n**************** LEAVING P4PFCODE ****************\n";
	
	return 0;
}
/*********************************************/
/*********** END fourth funciton *************/
/*********************************************/


/********************** C function *****************************/
// P4P + unknown focal length MEX helper
// given a set of 4x 2D<->3D correspondences, calculate action matrix
//
// by Martin Bujnak, (c)apr2008
//
//
// Please refer to the following paper, when using this code :
//
//      Bujnak, M., Kukelova, Z., and Pajdla, T. A general solution to the p4p
//      problem for camera with unknown focal length. CVPR 2008, Anchorage,
//      Alaska, USA, June 2008
//

double inline P4pf::dabs(double a)
{
	if (a > 0) return a;
	return -a;
}

void P4pf::GJ(double *A, int rcnt, int ccnt, double tol)
{
	int r = 0;      // row
	int c = 0;      // col
	int k;
	int l;
	int dstofs;
	int srcofs;
	int ofs = 0;
	int pofs = 0;
	double pivot_i;
	double b;
	
	// gj
	ofs = 0;
	pofs = 0;
	while (r < rcnt && c < ccnt) {
		
		// find pivot
		double apivot = 0;
		double pivot = 0;
		int pivot_r = -1;
		
		pofs = ofs;
		for (k = r; k < rcnt; k++) {
			
			// pivot selection criteria here !
			if (dabs(*(A+pofs)) > apivot) {
				
				pivot = *(A+pofs);
				apivot = dabs(pivot);
				pivot_r = k;
			}
			pofs += ccnt;
		}
		
		if (apivot < tol) {
			
			// empty col - shift to next col (or jump)
			c++;
			ofs++;
			
		} else {
			
			// process rows
			pivot_i = 1.0/pivot;
			
			// exchange pivot and selected rows
			// + divide row
			if (pivot_r == r) {
				
				srcofs = ofs;
				for (l = c; l < ccnt; l++) {
					
					*(A+srcofs) = *(A+srcofs)*pivot_i;
					srcofs++;
				}
				
			} else {
				
				srcofs = ofs;
				dstofs = ccnt*pivot_r+c;
				for (l = c; l < ccnt; l++) {
					
					b = *(A+srcofs);
					*(A+srcofs) = *(A+dstofs)*pivot_i;
					*(A+dstofs) = b;
					
					srcofs++;
					dstofs++;
				}
			}
			
			// zero bottom
			pofs = ofs + ccnt;
			for (k = r + 1; k < rcnt; k++) {
				
				if (dabs(*(A+pofs)) > tol) {
					
					// nonzero row
					b = *(A+pofs);
					dstofs = pofs + 1;
					srcofs = ofs + 1;
					for (l = c + 1; l < ccnt; l++) {
						
						*(A+dstofs) = (*(A+dstofs) - *(A+srcofs) * b);
						dstofs++;
						srcofs++;
					}
					*(A+pofs) = 0;
				}
				pofs += ccnt;
			}
			
			// zero top
			pofs = c;
			for (k = 0; k < r; k++) {
				
				if (dabs(*(A+pofs)) > tol) {
					
					// nonzero row
					b = *(A+pofs);
					dstofs = pofs + 1;
					srcofs = ofs + 1;
					for (l = c + 1; l < ccnt; l++) {
						
						*(A+dstofs) = (*(A+dstofs) - *(A+srcofs) * b);
						dstofs++;
						srcofs++;
					}
					*(A+pofs) = 0;
				}
				pofs += ccnt;
			}
			
			r++;
			c++;
			ofs += ccnt + 1;
		}
	}
}

//
// prepare polynomial coefficients

void P4pf::CalcCoefs(double const *src1, double const *src2, double const *src3, double const *src4, double const *src5, double *dst1)
{
	
	//	symbolic names.
	#define glab (src1[0])
	#define glac (src1[1])
	#define glad (src1[2])
	#define glbc (src1[3])
	#define glbd (src1[4])
	#define glcd (src1[5])
	#define a1 (src2[0])
	#define a2 (src2[1])
	#define b1 (src3[0])
	#define b2 (src3[1])
	#define c1 (src4[0])
	#define c2 (src4[1])
	#define d1 (src5[0])
	#define d2 (src5[1])
	
	double t1;
	double t11;
	double t12;
	double t13;
	double t14;
	double t16;
	double t2;
	double t24;
	double t27;
	double t28;
	double t3;
	double t32;
	double t33;
	double t34;
	double t35;
	double t37;
	double t39;
	double t4;
	double t41;
	double t42;
	double t43;
	double t46;
	double t5;
	double t51;
	double t53;
	double t56;
	double t59;
	double t60;
	double t67;
	double t84;
	double t9;
	
	//	consts
	
	
	
	// destination group 1
	
	t1 =  1/glad;
	t2 =  t1*glbc;
	t3 =  glab*t1;
	t4 =  glac*t1;
	t5 =  t2-t3-t4;
	t9 =  d2*d2;
	t11 =  t3*t9;
	t12 =  t4*t9;
	t13 =  d1*d1;
	t14 =  t4*t13;
	t16 =  t3*t13;
	t24 =  -a2*b2-b1*a1;
	t27 =  -c1*a1-c2*a2;
	t28 =  a1*t1;
	t32 =  t1*d1;
	t33 =  a1*glac*t32;
	t34 =  a2*d2;
	t35 =  t4*t34;
	t37 =  a1*glab*t32;
	t39 =  t3*t34;
	t41 =  a2*a2;
	t42 =  a1*a1;
	t43 =  t4*t41;
	t46 =  t42*glac*t1;
	t51 =  t42*glab*t1;
	t53 =  t42*t1;
	t56 =  t3*t41;
	t59 =  c1*c1;
	t60 =  c2*c2;
	t67 =  t1*glbd;
	t84 =  glcd*t1;
	
	
	// destination group 2
	
	
	
	// destination group 0
	
	
	
	// destination group 1
	
	dst1[0] = 1.0;
	dst1[1] = t5/2.0;
	dst1[2] = -1.0;
	dst1[3] = -1.0;
	dst1[4] = c2*b2+c1*b1;
	dst1[5] = -t5;
	dst1[6] = t2*t9/2.0-t11/2.0-t12/2.0-t14/2.0+t2*t13/2.0-t16/2.0;
	dst1[7] = 1.0-t4/2.0-t3/2.0+t2/2.0;
	dst1[8] = t24;
	dst1[9] = t27;
	dst1[10] = -t28*glbc*d1+t33+t35+t37-t2*t34+t39;
	dst1[11] = t41+t42-t43/2.0-t46/2.0+t2*t41/2.0-t51/2.0+t53*glbc/2.0-t56/2.0;
	dst1[12] = 1.0;
	dst1[13] = -t4;
	dst1[14] = -2.0;
	dst1[15] = t59+t60;
	dst1[16] = 2.0*t4;
	dst1[17] = -t14-t12;
	dst1[18] = -t4+1.0;
	dst1[19] = 2.0*t27;
	dst1[20] = 2.0*t33+2.0*t35;
	dst1[21] = -t43+t41+t42-t46;
	dst1[22] = 1.0;
	dst1[23] = t67/2.0-1.0/2.0-t3/2.0;
	dst1[24] = -1.0;
	dst1[25] = t3-t67;
	dst1[26] = d2*b2+b1*d1;
	dst1[27] = -t11/2.0-t16/2.0+t67*t9/2.0+t67*t13/2.0-t9/2.0-t13/2.0;
	dst1[28] = -t3/2.0+t67/2.0+1.0/2.0;
	dst1[29] = t24;
	dst1[30] = -t28*glbd*d1+t37+t39-t67*t34;
	dst1[31] = t67*t41/2.0+t53*glbd/2.0-t56/2.0-t51/2.0+t42/2.0+t41/2.0;
	dst1[32] = 1.0;
	dst1[33] = -t4/2.0+t84/2.0-1.0/2.0;
	dst1[34] = -1.0;
	dst1[35] = t4-t84;
	dst1[36] = c1*d1+c2*d2;
	dst1[37] = t84*t13/2.0-t12/2.0-t14/2.0-t13/2.0-t9/2.0+t84*t9/2.0;
	dst1[38] = -t4/2.0+1.0/2.0+t84/2.0;
	dst1[39] = t27;
	dst1[40] = t35+t33-t84*t34-glcd*a1*t32;
	dst1[41] = t42/2.0+t41/2.0-t43/2.0+glcd*t42*t1/2.0-t46/2.0+t84*t41/2.0;
	
	
	// destination group 2
	
	
	
	#undef glab
	#undef glac
	#undef glad
	#undef glbc
	#undef glbd
	#undef glcd
	#undef a1
	#undef a2
	#undef b1
	#undef b2
	#undef c1
	#undef c2
	#undef d1
	#undef d2
}

// [glab, glac, glad, glbc, glbd, glcd], [a1; a2], [b1; b2], [c1; c2], [d1;d2]
//  glXY - ||X-Y||^2 - quadratic distances between 3D points X and Y
//  a1 (a2) = x (resp y) measurement of the first 2D point
//  b1 (b2) = x (resp y) measurement of the second 2D point
//  c1 (c2) = x (resp y) measurement of the third 2D point
//  d1 (d2) = x (resp y) measurement of the fourth 2D point
//
// output *A - 10x10 action matrix

void P4pf::p4pfmex(double *glab, double *a1, double *b1, double *c1, double *d1, Eigen::Matrix<std::complex<double>, 10, 10> *A)
{
	// precalculate polynomial equations coefficients
	
	double M[6864];
	double coefs[42];
	
	P4pf::CalcCoefs(glab, a1, b1, c1, d1, coefs);
	
	//
	memset(M, 0, 6864*sizeof(double));
	
	M[64] = coefs[0];	M[403] = coefs[0];	M[486] = coefs[0];	M[572] = coefs[0];	M[1533] = coefs[0];	M[1616] = coefs[0];	M[1702] = coefs[0];	M[1787] = coefs[0];	M[1874] = coefs[0];	M[1960] = coefs[0];	M[3979] = coefs[0];	M[4063] = coefs[0];	M[4149] = coefs[0];	M[4234] = coefs[0];	M[4321] = coefs[0];	M[4407] = coefs[0];	M[4494] = coefs[0];	M[6161] = coefs[0];	M[6248] = coefs[0];
	M[71] = coefs[1];	M[411] = coefs[1];	M[496] = coefs[1];	M[582] = coefs[1];	M[1539] = coefs[1];	M[1626] = coefs[1];	M[1712] = coefs[1];	M[1798] = coefs[1];	M[1884] = coefs[1];	M[4071] = coefs[1];	M[4157] = coefs[1];	M[4244] = coefs[1];	M[4330] = coefs[1];	M[4416] = coefs[1];	M[4500] = coefs[1];	M[6165] = coefs[1];	M[6252] = coefs[1];
	M[75] = coefs[2];	M[419] = coefs[2];	M[504] = coefs[2];	M[590] = coefs[2];	M[1551] = coefs[2];	M[1635] = coefs[2];	M[1721] = coefs[2];	M[1806] = coefs[2];	M[1892] = coefs[2];	M[4001] = coefs[2];	M[4085] = coefs[2];	M[4171] = coefs[2];	M[4256] = coefs[2];	M[4342] = coefs[2];	M[4428] = coefs[2];	M[4512] = coefs[2];	M[6171] = coefs[2];	M[6255] = coefs[2];
	M[76] = coefs[3];	M[420] = coefs[3];	M[505] = coefs[3];	M[591] = coefs[3];	M[1552] = coefs[3];	M[1636] = coefs[3];	M[1722] = coefs[3];	M[1807] = coefs[3];	M[1893] = coefs[3];	M[4002] = coefs[3];	M[4086] = coefs[3];	M[4172] = coefs[3];	M[4257] = coefs[3];	M[4343] = coefs[3];	M[4429] = coefs[3];	M[4513] = coefs[3];	M[6172] = coefs[3];	M[6256] = coefs[3];
	M[77] = coefs[4];	M[421] = coefs[4];	M[506] = coefs[4];	M[592] = coefs[4];	M[1553] = coefs[4];	M[1637] = coefs[4];	M[1723] = coefs[4];	M[1808] = coefs[4];	M[1894] = coefs[4];	M[1980] = coefs[4];	M[4087] = coefs[4];	M[4173] = coefs[4];	M[4258] = coefs[4];	M[4344] = coefs[4];	M[4430] = coefs[4];	M[4514] = coefs[4];	M[6173] = coefs[4];	M[6257] = coefs[4];
	M[79] = coefs[5];	M[423] = coefs[5];	M[508] = coefs[5];	M[1555] = coefs[5];	M[1640] = coefs[5];	M[1726] = coefs[5];	M[1812] = coefs[5];	M[1898] = coefs[5];	M[4003] = coefs[5];	M[4090] = coefs[5];	M[4176] = coefs[5];	M[4262] = coefs[5];	M[4348] = coefs[5];	M[4517] = coefs[5];	M[6176] = coefs[5];	M[6260] = coefs[5];
	M[82] = coefs[6];	M[426] = coefs[6];	M[513] = coefs[6];	M[599] = coefs[6];	M[1645] = coefs[6];	M[1731] = coefs[6];	M[1818] = coefs[6];	M[1904] = coefs[6];	M[1990] = coefs[6];	M[4179] = coefs[6];	M[4354] = coefs[6];	M[4440] = coefs[6];	M[4524] = coefs[6];	M[6181] = coefs[6];	M[6266] = coefs[6];
	M[83] = coefs[7];	M[431] = coefs[7];	M[516] = coefs[7];	M[1567] = coefs[7];	M[1652] = coefs[7];	M[1825] = coefs[7];	M[1911] = coefs[7];	M[4019] = coefs[7];	M[4104] = coefs[7];	M[4190] = coefs[7];	M[4276] = coefs[7];	M[4362] = coefs[7];	M[6277] = coefs[7];
	M[84] = coefs[8];	M[432] = coefs[8];	M[517] = coefs[8];	M[603] = coefs[8];	M[1568] = coefs[8];	M[1653] = coefs[8];	M[1739] = coefs[8];	M[1826] = coefs[8];	M[1912] = coefs[8];	M[1998] = coefs[8];	M[4020] = coefs[8];	M[4105] = coefs[8];	M[4191] = coefs[8];	M[4277] = coefs[8];	M[4363] = coefs[8];	M[4449] = coefs[8];	M[4532] = coefs[8];	M[6195] = coefs[8];	M[6278] = coefs[8];
	M[85] = coefs[9];	M[433] = coefs[9];	M[518] = coefs[9];	M[604] = coefs[9];	M[1569] = coefs[9];	M[1654] = coefs[9];	M[1740] = coefs[9];	M[1913] = coefs[9];	M[1999] = coefs[9];	M[4021] = coefs[9];	M[4106] = coefs[9];	M[4192] = coefs[9];	M[4364] = coefs[9];	M[4450] = coefs[9];	M[4533] = coefs[9];	M[6196] = coefs[9];	M[6279] = coefs[9];
	M[86] = coefs[10];	M[434] = coefs[10];	M[521] = coefs[10];	M[607] = coefs[10];	M[1570] = coefs[10];	M[1657] = coefs[10];	M[1743] = coefs[10];	M[1830] = coefs[10];	M[1916] = coefs[10];	M[4109] = coefs[10];	M[4195] = coefs[10];	M[4282] = coefs[10];	M[4368] = coefs[10];	M[4454] = coefs[10];	M[4538] = coefs[10];	M[6200] = coefs[10];	M[6284] = coefs[10];
	M[87] = coefs[11];	M[438] = coefs[11];	M[525] = coefs[11];	M[611] = coefs[11];	M[1578] = coefs[11];	M[1665] = coefs[11];	M[1751] = coefs[11];	M[1838] = coefs[11];	M[1924] = coefs[11];	M[4034] = coefs[11];	M[4121] = coefs[11];	M[4207] = coefs[11];	M[4294] = coefs[11];	M[4380] = coefs[11];	M[4551] = coefs[11];	M[6214] = coefs[11];	M[6298] = coefs[11];
	M[153] = coefs[12];	M[668] = coefs[12];	M[750] = coefs[12];	M[837] = coefs[12];	M[2062] = coefs[12];	M[2145] = coefs[12];	M[2232] = coefs[12];	M[2319] = coefs[12];	M[2403] = coefs[12];	M[2490] = coefs[12];	M[2577] = coefs[12];	M[4596] = coefs[12];	M[4679] = coefs[12];	M[4766] = coefs[12];	M[4850] = coefs[12];	M[4937] = coefs[12];	M[5024] = coefs[12];	M[5110] = coefs[12];	M[6338] = coefs[12];	M[6424] = coefs[12];
	M[159] = coefs[13];	M[675] = coefs[13];	M[759] = coefs[13];	M[846] = coefs[13];	M[2067] = coefs[13];	M[2154] = coefs[13];	M[2241] = coefs[13];	M[2328] = coefs[13];	M[2413] = coefs[13];	M[2499] = coefs[13];	M[4686] = coefs[13];	M[4773] = coefs[13];	M[4859] = coefs[13];	M[4945] = coefs[13];	M[5032] = coefs[13];	M[5115] = coefs[13];	M[6341] = coefs[13];	M[6427] = coefs[13];
	M[164] = coefs[14];	M[684] = coefs[14];	M[768] = coefs[14];	M[855] = coefs[14];	M[2080] = coefs[14];	M[2164] = coefs[14];	M[2251] = coefs[14];	M[2338] = coefs[14];	M[2422] = coefs[14];	M[2508] = coefs[14];	M[4618] = coefs[14];	M[4701] = coefs[14];	M[4788] = coefs[14];	M[4872] = coefs[14];	M[4958] = coefs[14];	M[5045] = coefs[14];	M[5128] = coefs[14];	M[6348] = coefs[14];	M[6431] = coefs[14];
	M[166] = coefs[15];	M[686] = coefs[15];	M[770] = coefs[15];	M[857] = coefs[15];	M[2082] = coefs[15];	M[2253] = coefs[15];	M[2340] = coefs[15];	M[2424] = coefs[15];	M[2510] = coefs[15];	M[2597] = coefs[15];	M[4703] = coefs[15];	M[4790] = coefs[15];	M[4874] = coefs[15];	M[4960] = coefs[15];	M[5047] = coefs[15];	M[5130] = coefs[15];	M[6350] = coefs[15];	M[6433] = coefs[15];
	M[167] = coefs[16];	M[687] = coefs[16];	M[771] = coefs[16];	M[2083] = coefs[16];	M[2168] = coefs[16];	M[2255] = coefs[16];	M[2342] = coefs[16];	M[2427] = coefs[16];	M[2513] = coefs[16];	M[4619] = coefs[16];	M[4705] = coefs[16];	M[4792] = coefs[16];	M[4877] = coefs[16];	M[4963] = coefs[16];	M[5132] = coefs[16];	M[6352] = coefs[16];	M[6435] = coefs[16];
	M[170] = coefs[17];	M[690] = coefs[17];	M[776] = coefs[17];	M[863] = coefs[17];	M[2173] = coefs[17];	M[2260] = coefs[17];	M[2347] = coefs[17];	M[2433] = coefs[17];	M[2519] = coefs[17];	M[2606] = coefs[17];	M[4795] = coefs[17];	M[4969] = coefs[17];	M[5056] = coefs[17];	M[5139] = coefs[17];	M[6357] = coefs[17];	M[6441] = coefs[17];
	M[171] = coefs[18];	M[695] = coefs[18];	M[779] = coefs[18];	M[2095] = coefs[18];	M[2180] = coefs[18];	M[2267] = coefs[18];	M[2440] = coefs[18];	M[2526] = coefs[18];	M[4635] = coefs[18];	M[4719] = coefs[18];	M[4806] = coefs[18];	M[4891] = coefs[18];	M[4977] = coefs[18];	M[6452] = coefs[18];
	M[173] = coefs[19];	M[697] = coefs[19];	M[781] = coefs[19];	M[868] = coefs[19];	M[2097] = coefs[19];	M[2182] = coefs[19];	M[2269] = coefs[19];	M[2356] = coefs[19];	M[2442] = coefs[19];	M[2528] = coefs[19];	M[2615] = coefs[19];	M[4637] = coefs[19];	M[4721] = coefs[19];	M[4808] = coefs[19];	M[4893] = coefs[19];	M[4979] = coefs[19];	M[5066] = coefs[19];	M[5148] = coefs[19];	M[6372] = coefs[19];	M[6454] = coefs[19];
	M[174] = coefs[20];	M[698] = coefs[20];	M[784] = coefs[20];	M[871] = coefs[20];	M[2098] = coefs[20];	M[2185] = coefs[20];	M[2272] = coefs[20];	M[2359] = coefs[20];	M[2445] = coefs[20];	M[2531] = coefs[20];	M[4724] = coefs[20];	M[4811] = coefs[20];	M[4897] = coefs[20];	M[4983] = coefs[20];	M[5070] = coefs[20];	M[5153] = coefs[20];	M[6376] = coefs[20];	M[6459] = coefs[20];
	M[175] = coefs[21];	M[702] = coefs[21];	M[788] = coefs[21];	M[875] = coefs[21];	M[2106] = coefs[21];	M[2193] = coefs[21];	M[2280] = coefs[21];	M[2367] = coefs[21];	M[2453] = coefs[21];	M[2539] = coefs[21];	M[4650] = coefs[21];	M[4736] = coefs[21];	M[4823] = coefs[21];	M[4909] = coefs[21];	M[4995] = coefs[21];	M[5166] = coefs[21];	M[6390] = coefs[21];	M[6473] = coefs[21];
	M[243] = coefs[22];	M[935] = coefs[22];	M[1019] = coefs[22];	M[1105] = coefs[22];	M[2681] = coefs[22];	M[2765] = coefs[22];	M[2851] = coefs[22];	M[2936] = coefs[22];	M[3022] = coefs[22];	M[3108] = coefs[22];	M[5209] = coefs[22];	M[5293] = coefs[22];	M[5379] = coefs[22];	M[5463] = coefs[22];	M[6515] = coefs[22];	M[6601] = coefs[22];
	M[247] = coefs[23];	M[939] = coefs[23];	M[1024] = coefs[23];	M[1110] = coefs[23];	M[2683] = coefs[23];	M[2770] = coefs[23];	M[2856] = coefs[23];	M[2942] = coefs[23];	M[3028] = coefs[23];	M[5213] = coefs[23];	M[5298] = coefs[23];	M[5384] = coefs[23];	M[5468] = coefs[23];	M[6517] = coefs[23];	M[6604] = coefs[23];
	M[251] = coefs[24];	M[947] = coefs[24];	M[1032] = coefs[24];	M[1118] = coefs[24];	M[2695] = coefs[24];	M[2779] = coefs[24];	M[2865] = coefs[24];	M[2950] = coefs[24];	M[3036] = coefs[24];	M[5227] = coefs[24];	M[5310] = coefs[24];	M[5396] = coefs[24];	M[5480] = coefs[24];	M[6523] = coefs[24];	M[6607] = coefs[24];
	M[255] = coefs[25];	M[951] = coefs[25];	M[1036] = coefs[25];	M[2699] = coefs[25];	M[2784] = coefs[25];	M[2870] = coefs[25];	M[2956] = coefs[25];	M[3042] = coefs[25];	M[5232] = coefs[25];	M[5316] = coefs[25];	M[5485] = coefs[25];	M[6528] = coefs[25];	M[6612] = coefs[25];
	M[256] = coefs[26];	M[952] = coefs[26];	M[1037] = coefs[26];	M[1123] = coefs[26];	M[2700] = coefs[26];	M[2785] = coefs[26];	M[2871] = coefs[26];	M[2957] = coefs[26];	M[3043] = coefs[26];	M[3129] = coefs[26];	M[5233] = coefs[26];	M[5317] = coefs[26];	M[5403] = coefs[26];	M[5486] = coefs[26];	M[6529] = coefs[26];	M[6613] = coefs[26];
	M[258] = coefs[27];	M[954] = coefs[27];	M[1041] = coefs[27];	M[1127] = coefs[27];	M[2789] = coefs[27];	M[2875] = coefs[27];	M[2962] = coefs[27];	M[3048] = coefs[27];	M[3134] = coefs[27];	M[5235] = coefs[27];	M[5322] = coefs[27];	M[5408] = coefs[27];	M[5492] = coefs[27];	M[6533] = coefs[27];	M[6618] = coefs[27];
	M[259] = coefs[28];	M[959] = coefs[28];	M[1044] = coefs[28];	M[2711] = coefs[28];	M[2796] = coefs[28];	M[2969] = coefs[28];	M[3055] = coefs[28];	M[5246] = coefs[28];	M[5330] = coefs[28];	M[6629] = coefs[28];
	M[260] = coefs[29];	M[960] = coefs[29];	M[1045] = coefs[29];	M[1131] = coefs[29];	M[2712] = coefs[29];	M[2797] = coefs[29];	M[2883] = coefs[29];	M[2970] = coefs[29];	M[3056] = coefs[29];	M[3142] = coefs[29];	M[5247] = coefs[29];	M[5331] = coefs[29];	M[5417] = coefs[29];	M[5500] = coefs[29];	M[6547] = coefs[29];	M[6630] = coefs[29];
	M[262] = coefs[30];	M[962] = coefs[30];	M[1049] = coefs[30];	M[1135] = coefs[30];	M[2714] = coefs[30];	M[2801] = coefs[30];	M[2887] = coefs[30];	M[2974] = coefs[30];	M[3060] = coefs[30];	M[5251] = coefs[30];	M[5336] = coefs[30];	M[5422] = coefs[30];	M[5506] = coefs[30];	M[6552] = coefs[30];	M[6636] = coefs[30];
	M[263] = coefs[31];	M[966] = coefs[31];	M[1053] = coefs[31];	M[1139] = coefs[31];	M[2722] = coefs[31];	M[2809] = coefs[31];	M[2895] = coefs[31];	M[2982] = coefs[31];	M[3068] = coefs[31];	M[5263] = coefs[31];	M[5348] = coefs[31];	M[5519] = coefs[31];	M[6566] = coefs[31];	M[6650] = coefs[31];
	M[332] = coefs[32];	M[1200] = coefs[32];	M[1284] = coefs[32];	M[1371] = coefs[32];	M[1458] = coefs[32];	M[3210] = coefs[32];	M[3294] = coefs[32];	M[3381] = coefs[32];	M[3468] = coefs[32];	M[3553] = coefs[32];	M[3640] = coefs[32];	M[3727] = coefs[32];	M[3814] = coefs[32];	M[3901] = coefs[32];	M[5564] = coefs[32];	M[5651] = coefs[32];	M[5738] = coefs[32];	M[5822] = coefs[32];	M[5908] = coefs[32];	M[5994] = coefs[32];	M[6080] = coefs[32];	M[6692] = coefs[32];	M[6778] = coefs[32];
	M[335] = coefs[33];	M[1203] = coefs[33];	M[1288] = coefs[33];	M[1375] = coefs[33];	M[1462] = coefs[33];	M[3211] = coefs[33];	M[3298] = coefs[33];	M[3385] = coefs[33];	M[3472] = coefs[33];	M[3558] = coefs[33];	M[3645] = coefs[33];	M[3732] = coefs[33];	M[3819] = coefs[33];	M[5567] = coefs[33];	M[5654] = coefs[33];	M[5741] = coefs[33];	M[5826] = coefs[33];	M[5912] = coefs[33];	M[5999] = coefs[33];	M[6084] = coefs[33];	M[6693] = coefs[33];	M[6780] = coefs[33];
	M[340] = coefs[34];	M[1212] = coefs[34];	M[1297] = coefs[34];	M[1384] = coefs[34];	M[1471] = coefs[34];	M[3224] = coefs[34];	M[3308] = coefs[34];	M[3395] = coefs[34];	M[3482] = coefs[34];	M[3567] = coefs[34];	M[3654] = coefs[34];	M[3741] = coefs[34];	M[3828] = coefs[34];	M[5582] = coefs[34];	M[5669] = coefs[34];	M[5756] = coefs[34];	M[5839] = coefs[34];	M[5925] = coefs[34];	M[6011] = coefs[34];	M[6097] = coefs[34];	M[6700] = coefs[34];	M[6784] = coefs[34];
	M[343] = coefs[35];	M[1215] = coefs[35];	M[1300] = coefs[35];	M[1387] = coefs[35];	M[3227] = coefs[35];	M[3312] = coefs[35];	M[3399] = coefs[35];	M[3486] = coefs[35];	M[3572] = coefs[35];	M[3659] = coefs[35];	M[3746] = coefs[35];	M[3833] = coefs[35];	M[5586] = coefs[35];	M[5673] = coefs[35];	M[5760] = coefs[35];	M[5844] = coefs[35];	M[6016] = coefs[35];	M[6101] = coefs[35];	M[6704] = coefs[35];	M[6788] = coefs[35];
	M[345] = coefs[36];	M[1217] = coefs[36];	M[1302] = coefs[36];	M[1389] = coefs[36];	M[1476] = coefs[36];	M[3229] = coefs[36];	M[3314] = coefs[36];	M[3401] = coefs[36];	M[3488] = coefs[36];	M[3661] = coefs[36];	M[3748] = coefs[36];	M[3835] = coefs[36];	M[3922] = coefs[36];	M[5762] = coefs[36];	M[5846] = coefs[36];	M[5932] = coefs[36];	M[6018] = coefs[36];	M[6103] = coefs[36];	M[6706] = coefs[36];	M[6790] = coefs[36];
	M[346] = coefs[37];	M[1218] = coefs[37];	M[1305] = coefs[37];	M[1392] = coefs[37];	M[1479] = coefs[37];	M[3317] = coefs[37];	M[3404] = coefs[37];	M[3491] = coefs[37];	M[3578] = coefs[37];	M[3665] = coefs[37];	M[3752] = coefs[37];	M[3839] = coefs[37];	M[3926] = coefs[37];	M[5763] = coefs[37];	M[5850] = coefs[37];	M[5936] = coefs[37];	M[6023] = coefs[37];	M[6108] = coefs[37];	M[6709] = coefs[37];	M[6794] = coefs[37];
	M[347] = coefs[38];	M[1223] = coefs[38];	M[1308] = coefs[38];	M[1395] = coefs[38];	M[3239] = coefs[38];	M[3324] = coefs[38];	M[3411] = coefs[38];	M[3585] = coefs[38];	M[3672] = coefs[38];	M[3759] = coefs[38];	M[3846] = coefs[38];	M[5600] = coefs[38];	M[5687] = coefs[38];	M[5774] = coefs[38];	M[5858] = coefs[38];	M[6030] = coefs[38];	M[6805] = coefs[38];
	M[349] = coefs[39];	M[1225] = coefs[39];	M[1310] = coefs[39];	M[1397] = coefs[39];	M[1484] = coefs[39];	M[3241] = coefs[39];	M[3326] = coefs[39];	M[3413] = coefs[39];	M[3500] = coefs[39];	M[3674] = coefs[39];	M[3761] = coefs[39];	M[3848] = coefs[39];	M[3935] = coefs[39];	M[5602] = coefs[39];	M[5689] = coefs[39];	M[5776] = coefs[39];	M[5860] = coefs[39];	M[5946] = coefs[39];	M[6032] = coefs[39];	M[6117] = coefs[39];	M[6724] = coefs[39];	M[6807] = coefs[39];
	M[350] = coefs[40];	M[1226] = coefs[40];	M[1313] = coefs[40];	M[1400] = coefs[40];	M[1487] = coefs[40];	M[3242] = coefs[40];	M[3329] = coefs[40];	M[3416] = coefs[40];	M[3503] = coefs[40];	M[3590] = coefs[40];	M[3677] = coefs[40];	M[3764] = coefs[40];	M[3851] = coefs[40];	M[5605] = coefs[40];	M[5692] = coefs[40];	M[5779] = coefs[40];	M[5864] = coefs[40];	M[5950] = coefs[40];	M[6037] = coefs[40];	M[6122] = coefs[40];	M[6728] = coefs[40];	M[6812] = coefs[40];
	M[351] = coefs[41];	M[1230] = coefs[41];	M[1317] = coefs[41];	M[1404] = coefs[41];	M[1491] = coefs[41];	M[3250] = coefs[41];	M[3337] = coefs[41];	M[3424] = coefs[41];	M[3511] = coefs[41];	M[3598] = coefs[41];	M[3685] = coefs[41];	M[3772] = coefs[41];	M[3859] = coefs[41];	M[5617] = coefs[41];	M[5704] = coefs[41];	M[5791] = coefs[41];	M[5876] = coefs[41];	M[6050] = coefs[41];	M[6135] = coefs[41];	M[6742] = coefs[41];	M[6826] = coefs[41];
	
	
	// GJ elimination
	P4pf::GJ(M, 78, 88, 2.2204e-11);
	
	
	// action matrix
	memset(A, 0, sizeof(double)*100);
	
	// Modified these in order to fit in the C++ code.
	(*A)(1) = 1;
	(*A)(15) = 1;
	(*A)(26) = 1;
	(*A)(37) = 1;
	(*A)(48) = 1;
	(*A)(50) = -M[6599];	(*A)(51) = -M[6598];	(*A)(52) = -M[6597];	(*A)(53) = -M[6596];	(*A)(54) = -M[6595];	(*A)(55) = -M[6594];	(*A)(56) = -M[6593];	(*A)(57) = -M[6592];	(*A)(58) = -M[6591];	(*A)(59) = -M[6590];
	(*A)(60) = -M[6511];	(*A)(61) = -M[6510];	(*A)(62) = -M[6509];	(*A)(63) = -M[6508];	(*A)(64) = -M[6507];	(*A)(65) = -M[6506];	(*A)(66) = -M[6505];	(*A)(67) = -M[6504];	(*A)(68) = -M[6503];	(*A)(69) = -M[6502];
	(*A)(70) = -M[6423];	(*A)(71) = -M[6422];	(*A)(72) = -M[6421];	(*A)(73) = -M[6420];	(*A)(74) = -M[6419];	(*A)(75) = -M[6418];	(*A)(76) = -M[6417];	(*A)(77) = -M[6416];	(*A)(78) = -M[6415];	(*A)(79) = -M[6414];
	(*A)(80) = -M[6335];	(*A)(81) = -M[6334];	(*A)(82) = -M[6333];	(*A)(83) = -M[6332];	(*A)(84) = -M[6331];	(*A)(85) = -M[6330];	(*A)(86) = -M[6329];	(*A)(87) = -M[6328];	(*A)(88) = -M[6327];	(*A)(89) = -M[6326];
	(*A)(90) = -M[6247];	(*A)(91) = -M[6246];	(*A)(92) = -M[6245];	(*A)(93) = -M[6244];	(*A)(94) = -M[6243];	(*A)(95) = -M[6242];	(*A)(96) = -M[6241];	(*A)(97) = -M[6240];	(*A)(98) = -M[6239];	(*A)(99) = -M[6238];
	
}


